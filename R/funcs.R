#' Query NOAA VDatum API for MLLW elevation in NAVD88
#'
#' Converts 0 m MLLW to NAVD88 at a single point using the NOAA VDatum web API.
#'
#' @param lon Numeric. Longitude in decimal degrees (NAD83).
#' @param lat Numeric. Latitude in decimal degrees (NAD83).
#' @param unit Character. Vertical unit for input/output, default \code{"m"} for meters.
#'
#' @return Numeric. NAVD88 elevation of MLLW in the requested unit, or \code{NA} if
#'   the point is outside VDatum coverage or the request fails.

query_vdatum <- function(lon, lat, unit = "m") {
  resp <- tryCatch(
    {
      httr2::request("https://vdatum.noaa.gov/vdatumweb/api/convert") |>
        httr2::req_url_query(
          s_x = round(lon, 6),
          s_y = round(lat, 6),
          s_h_frame = "NAD83_2011",
          s_v_frame = "MLLW",
          s_v_unit = unit,
          s_v_geoid = "geoid18",
          t_h_frame = "NAD83_2011",
          t_v_frame = "NAVD88",
          t_v_unit = unit,
          t_v_geoid = "geoid18",
          s_v_value = 0
        ) |>
        httr2::req_timeout(15) |>
        httr2::req_retry(max_tries = 3, backoff = ~2) |>
        httr2::req_perform()
    },
    error = function(e) {
      message("Request error at (", lon, ", ", lat, "): ", e$message)
      NULL
    }
  )

  if (is.null(resp)) return(NA_real_)

  result <- tryCatch(httr2::resp_body_json(resp), error = function(e) NULL)
  if (is.null(result)) return(NA_real_)

  val <- suppressWarnings(as.numeric(result$t_z))

  if (!is.na(val) && val != -999999) return(val)

  NA_real_
}

#' Build a regular VDatum sample grid clipped to the study area
#'
#' Creates a regular lon/lat grid over the bounding box and clips it to the
#' county boundaries, for use with \code{batch_query_vdatum()}.
#'
#' @param bbox_4326 A named numeric vector from \code{sf::st_bbox()} in EPSG:4326.
#' @param counties_4326 An \code{sf} polygon of the study area counties in EPSG:4326.
#' @param grid_spacing Numeric. Grid spacing in decimal degrees. Default \code{0.05}
#'   (~5 km at this latitude).
#'
#' @return A tibble with columns \code{lon} and \code{lat} for points inside the study area.

build_vdatum_grid <- function(bbox_4326, counties_4326, grid_spacing = 0.05) {
  grid_pts <- expand.grid(
    lon = seq(bbox_4326["xmin"], bbox_4326["xmax"], by = grid_spacing),
    lat = seq(bbox_4326["ymin"], bbox_4326["ymax"], by = grid_spacing)
  ) |>
    dplyr::as_tibble()

  grid_sf <- sf::st_as_sf(grid_pts, coords = c("lon", "lat"), crs = 4326)
  in_area <- sf::st_intersects(grid_sf, sf::st_union(counties_4326), sparse = FALSE)[, 1]
  grid_pts[in_area, ]
}

#' Batch query VDatum API over a sample grid
#'
#' Calls \code{query_vdatum()} for every row in \code{grid_pts} with a 0.5-second
#' delay between requests (~2 req/sec) to respect NOAA rate limits.
#'
#' @param grid_pts A tibble with columns \code{lon} and \code{lat} in decimal degrees
#'   (NAD83), as returned by \code{build_vdatum_grid()}.
#'
#' @return The input tibble with an additional column \code{mllw_navd88_m} containing
#'   the NAVD88 elevation of MLLW in meters. Points outside coverage are \code{NA}.

batch_query_vdatum <- function(grid_pts) {
  offsets <- vector("numeric", nrow(grid_pts))

  for (i in seq_len(nrow(grid_pts))) {
    offsets[i] <- query_vdatum(grid_pts$lon[i], grid_pts$lat[i])

    if (i %% 10 == 0 || i == nrow(grid_pts)) {
      cat(sprintf("  %d / %d complete (%.0f%%)\n", i, nrow(grid_pts), 100 * i / nrow(grid_pts)))
    }

    Sys.sleep(0.5) # ~2 req/sec
  }

  grid_pts$mllw_navd88_m <- offsets
  grid_pts
}

#' Interpolate a spatially varying MLLW offset surface
#'
#' Fits a thin-plate spline to the VDatum sample points and predicts onto a
#' fine grid, then projects to EPSG:3087 and resamples to exactly match the
#' DEM grid.
#'
#' @param grid_pts A tibble with columns \code{lon}, \code{lat}, and
#'   \code{mllw_navd88_m}, as returned by \code{batch_query_vdatum()}. Rows with
#'   \code{NA} in \code{mllw_navd88_m} are dropped before fitting.
#' @param bbox_4326 A named numeric vector from \code{sf::st_bbox()} in EPSG:4326,
#'   used to define the prediction grid extent.
#' @param counties_3087 An \code{sf} polygon of the study area counties in EPSG:3087,
#'   used to mask the interpolated surface.
#' @param dem A \code{SpatRaster} in EPSG:3087. The output surface will be resampled
#'   to match this raster's grid exactly.
#'
#' @return A \code{SpatRaster} in EPSG:3087 with the same extent, resolution, and
#'   CRS as \code{dem}, containing MLLW elevation in meters NAVD88.

build_mllw_surface <- function(grid_pts, bbox_4326, counties_3087, dem) {
  pts_clean <- tidyr::drop_na(grid_pts, mllw_navd88_m)

  cat(sprintf("Interpolating from %d valid points...\n", nrow(pts_clean)))
  cat(sprintf(
    "MLLW offset range: %.3f to %.3f m NAVD88\n",
    min(pts_clean$mllw_navd88_m),
    max(pts_clean$mllw_navd88_m)
  ))

  tps_model <- fields::Tps(
    x = as.matrix(pts_clean[, c("lon", "lat")]),
    Y = pts_clean$mllw_navd88_m
  )

  interp_res <- 0.005 # ~500 m at this latitude
  pred_grid <- expand.grid(
    lon = seq(bbox_4326["xmin"], bbox_4326["xmax"], by = interp_res),
    lat = seq(bbox_4326["ymin"], bbox_4326["ymax"], by = interp_res)
  )
  pred_grid$mllw_offset_m <- predict(tps_model, as.matrix(pred_grid))

  offset_rast_4326 <- terra::rast(pred_grid, type = "xyz", crs = "EPSG:4326")
  offset_rast_3087 <- terra::project(offset_rast_4326, "EPSG:3087", method = "bilinear")
  offset_rast_3087 <- terra::mask(offset_rast_3087, terra::vect(counties_3087))

  mllw_surface <- terra::resample(offset_rast_3087, dem, method = "bilinear")

  stopifnot(terra::compareGeom(dem, mllw_surface, stopOnError = FALSE))

  mllw_surface
}

#' Delineate the coastal stratum polygon
#'
#' Identifies raster cells between MLLW and an upper elevation bound, then
#' converts to a dissolved \code{sf} polygon with stratum attributes.
#'
#' @param dem A \code{SpatRaster} in EPSG:3087 with elevation in meters NAVD88.
#' @param mllw_surface A \code{SpatRaster} aligned to \code{dem} containing the
#'   NAVD88 elevation of MLLW in meters, as returned by \code{build_mllw_surface()}.
#' @param upper_bound_m Numeric. Upper elevation limit of the stratum in meters
#'   NAVD88. Default is \code{5 * 0.3048} (5 ft = 1.524 m).
#'
#' @return An \code{sf} polygon in EPSG:3087 with columns \code{stratum},
#'   \code{lower_datum}, \code{upper_elev_ft}, \code{upper_elev_m}, and \code{crs}.

make_coastal_stratum <- function(dem, mllw_surface, upper_bound_m = 5 * 0.3048) {
  cat("Creating coastal stratum layer...\n")

  stratum_rast <- dem >= mllw_surface & dem <= upper_bound_m
  stratum_rast[!stratum_rast] <- NA

  cat("Converting to vector...\n")
  terra::as.polygons(stratum_rast, dissolve = TRUE) |>
    sf::st_as_sf() |>
    sf::st_make_valid() |>
    dplyr::mutate(
      stratum = "coastal",
      lower_datum = "MLLW (spatially varying, NAVD88)",
      upper_elev_ft = 5,
      upper_elev_m = upper_bound_m,
      crs = "EPSG:3087"
    )
}

#' Build a classified soils layer for a multi-county study area
#'
#' Fetches SSURGO map unit polygons county by county using
#' \code{fetch_soils_tiled()}, classifies each unit as Xeric, Mesic, or Hydric
#' following Ries and Scheda (2014) via \code{classify_soils_muname()}, clips
#' to the study area boundary, and dissolves by gridcode.
#'
#' @param counties An \code{sf} polygon with one row per county. Each row is
#'   passed individually to \code{fetch_soils_tiled()}. If a \code{county}
#'   column is present its values are used in progress messages.
#' @param crs Integer EPSG code for the output CRS. Default \code{3087}.
#' @param verbose Logical. If \code{TRUE} (default), prints per-county fetch
#'   progress and a summary of unclassified map units after classification.
#'   Water bodies will always appear in the unclassified list and can be
#'   ignored; any other entries may need adding to the series lists inside
#'   \code{classify_soils_muname()}.
#'
#' @return An \code{sf} object in \code{crs} with columns \code{gridcode}
#'   (100L = Xeric, 200L = Mesic, 300L = Hydric) and \code{Descrip}, dissolved
#'   and clipped to the union of \code{counties}.

build_soils_layer <- function(counties, crs = 3087L, verbose = TRUE) {

  # --- step 1: fetch polygons county by county -----------------------------

  mupolygon_list <- vector('list', nrow(counties))

  for (i in seq_len(nrow(counties))) {
    if (verbose) {
      label <- if ('county' %in% names(counties)) counties$county[i] else i
      cat('Fetching soil polygons for', label, '...\n')
    }
    mupolygon_list[[i]] <- fetch_soils_tiled(counties[i, ])
  }

  mupolygon <- dplyr::bind_rows(mupolygon_list)

  # --- steps 2-3: query munames and classify -------------------------------

  comp_classified <- classify_soils_muname(unique(mupolygon$mukey))

  if (verbose) {
    n_unclass <- sum(is.na(comp_classified$gridcode))
    cat('Unclassified mukeys:', n_unclass, '\n')
    if (n_unclass > 0) {
      comp_classified |>
        dplyr::filter(is.na(gridcode)) |>
        dplyr::count(muname, sort = TRUE) |>
        print()
    }
  }

  # --- step 4: join, clip to study area, dissolve --------------------------

  area_union <- sf::st_union(counties)

  mupolygon |>
    dplyr::left_join(
      comp_classified |> dplyr::select(mukey, gridcode, Descrip),
      by = 'mukey'
    ) |>
    dplyr::filter(!is.na(gridcode)) |>
    sf::st_transform(crs) |>
    sf::st_make_valid() |>
    sf::st_intersection(area_union) |>
    dplyr::group_by(gridcode, Descrip) |>
    dplyr::summarise(.groups = 'drop') |>
    dplyr::arrange(gridcode)
}

#' Classify SSURGO map units as Xeric, Mesic, or Hydric
#'
#' Queries map unit names from NRCS Soil Data Access (SDA) for a vector of
#' mukeys and classifies each unit using soil series name matching, following
#' the methodology of Ries and Scheda (2014). Xeric units are checked first so
#' that complex map unit names (e.g. "Tavares-Myakka complex") resolve to the
#' leading series. Hydric qualifiers in the map unit name (muck, depressional,
#' mostly wetland) override series-level matches to handle series like Pompano
#' and Basinger that appear in both Mesic and Hydric categories.
#'
#' @param mukeys Character or integer vector of SSURGO map unit keys.
#' @param chunk_size Integer. Number of mukeys per SDA query. Default \code{500}.
#'
#' @return A \code{data.frame} with columns \code{mukey}, \code{muname},
#'   \code{gridcode} (100L = Xeric, 200L = Mesic, 300L = Hydric, \code{NA} =
#'   unclassified), and \code{Descrip}. Water bodies and unrecognised map units
#'   will have \code{NA} gridcode and can be dropped downstream.

classify_soils_muname <- function(mukeys, chunk_size = 500L) {

  # --- series lookup lists (Ries & Scheda 2014) ----------------------------

  series_hydric <- c(
    'Sellers', 'Zephyr', 'Lacoochee', 'Okeelanta', 'Terra Ceia',
    'Weekiwachee', 'Samsula', 'Homosassa', 'Tobiska', 'Bessie', 'Delray',
    'Floridana', 'Kesson', 'Wulfert', 'Manatee', 'Okeechobee', 'Palmetto',
    'Broward', 'Matlacha', 'Haplaquents', 'Hydraquents', 'Holopaw', 'Nittaw',
    'Copeland', 'Waccasassa'
  )

  series_mesic <- c(
    'Wauchula', 'Pomona', 'Pineda', 'Felda', 'Myakka', 'Ora', 'Vero',
    'Immokalee', 'Smyrna', 'Basinger', 'Anclote', 'Pompano', 'EauGallie',
    'Eau Gallie', 'Chobee', 'Paisley', 'Arredondo', 'Cassia', 'Blichton',
    'Flemington', 'Placid', 'Aripeka', 'Wabasso', 'Ona', 'Pinellas',
    'St. Johns', 'Winder', 'Bradenton', 'Canova', 'Malabar', 'Seffner',
    'Parkwood', 'Braden', 'Eaton', 'Farmton', 'Kanapaha', 'Lynne',
    'Nobleton', 'Oldsmar', 'Waveland', 'Brynwood', 'Jumper', 'Lutterloh',
    'Mabel', 'Masaryk', 'Pople', 'Punta', 'Tarrytown', 'Wekiva'
  )

  series_xeric <- c(
    'Tavares', 'Sparr', 'Adamsville', 'Astatula', 'Chandler', 'Electra',
    'Paola', 'Narcoossee', 'Lake', 'Pomello', 'Kendrick', 'Lochloosa',
    'Newnan', 'Gainesville', 'Micanopy', 'Millhopper', 'Orlando', 'Zofo',
    'Zolfo', 'Candler', 'Nobleston', 'Beaches', 'Pits', 'Urban land',
    'Quartzipsamments', 'Arents', 'Dumps', 'Canaveral', 'Orsino',
    'Palm Beach', 'St. Augustine', 'Apopka', 'Archbold', 'Duette',
    'Fort Meade', 'Gypsum land', 'Jonesville', 'Neilhurst', 'Udorthents',
    'Slickens', 'Citronelle', 'Florahome', 'Ft. Green', 'Jonathan',
    'Redlevel', 'Shadeville', 'St. Lucie', 'Sumterville', 'Williston'
  )

  pat_hydric <- paste0('\\b(', paste(series_hydric, collapse = '|'), ')\\b')
  pat_mesic  <- paste0('\\b(', paste(series_mesic,  collapse = '|'), ')\\b')
  pat_xeric  <- paste0('\\b(', paste(series_xeric,  collapse = '|'), ')\\b')

  # --- query munames from SDA in chunks ------------------------------------

  chunks <- split(mukeys, ceiling(seq_along(mukeys) / chunk_size))

  comp_data <- lapply(chunks, function(keys) {
    soilDB::SDA_query(sprintf(
      "SELECT mu.mukey, mu.muname
         FROM mapunit mu
        WHERE mu.mukey IN (%s)",
      paste(keys, collapse = ',')
    ))
  }) |>
    dplyr::bind_rows()

  # --- classify by map unit name -------------------------------------------

  comp_data |>
    dplyr::mutate(
      gridcode = dplyr::case_when(
        stringr::str_detect(muname, stringr::regex(pat_xeric,  ignore_case = TRUE)) ~ 100L,
        stringr::str_detect(muname, stringr::regex('muck|depressional|mostly wetland', ignore_case = TRUE)) ~ 300L,
        stringr::str_detect(muname, stringr::regex(pat_hydric, ignore_case = TRUE)) ~ 300L,
        stringr::str_detect(muname, stringr::regex(pat_mesic,  ignore_case = TRUE)) ~ 200L,
        TRUE ~ NA_integer_
      ),
      Descrip = dplyr::case_when(
        gridcode == 100L ~ 'Xeric',
        gridcode == 200L ~ 'Mesic',
        gridcode == 300L ~ 'Hydric',
        TRUE             ~ 'Unclassified'
      )
    )
}

#' Fetch SSURGO map unit polygons using tiled SDA queries
#'
#' Subdivides the input geometry into a regular grid of tiles and issues one
#' \code{SDA_spatialQuery()} per tile, then binds all results. This avoids the
#' 90-second SDA timeout that occurs when querying large county geometries in
#' a single request.
#'
#' @param geom An \code{sf} polygon in any CRS; internally transformed to EPSG:4326.
#' @param tile_size Numeric. Tile size in decimal degrees. Default \code{0.1} (~10 km).
#' @param pause Numeric. Seconds to pause between tile requests. Default \code{0.5}.
#'
#' @return An \code{sf} object with columns \code{mukey}, \code{area_ac}, and geometry.
#'   Rows are not deduplicated — polygons spanning tile boundaries appear as multiple
#'   fragments, which are reconciled when dissolving by soil type downstream.

fetch_soils_tiled <- function(geom, tile_size = 0.1, pause = 0.5) {

  geom_4326 <- sf::st_transform(geom, 4326)

  tiles <- sf::st_make_grid(geom_4326, cellsize = tile_size) |>
    sf::st_as_sf() |>
    sf::st_filter(geom_4326)

  cat(sprintf('  Querying %d tiles (%.2f deg each)...\n', nrow(tiles), tile_size))

  results <- vector('list', nrow(tiles))

  for (j in seq_len(nrow(tiles))) {
    results[[j]] <- tryCatch(
      soilDB::SDA_spatialQuery(
        geom             = tiles[j, ],
        what             = 'mupolygon',
        db               = 'SSURGO',
        geomIntersection = TRUE
      ),
      error = function(e) {
        message(sprintf('  Tile %d/%d failed: %s', j, nrow(tiles), conditionMessage(e)))
        NULL
      }
    )

    Sys.sleep(pause)
  }

  dplyr::bind_rows(Filter(Negate(is.null), results))

}

#' Download a single WQP endpoint chunk to a cached zip file
#'
#' Issues a GET request to either the WQP Station or Result search endpoint for
#' one spatial unit (bounding box or county code), saves the compressed response
#' to \code{cache_dir}, and returns the zip file path. Re-running with the same
#' arguments skips the download and returns the cached path immediately.
#'
#' @param endpoint Character. WQP endpoint name: \code{"Station"} or \code{"Result"}.
#' @param label Character. Short label used in the cache filename and progress
#'   messages (e.g. a county name or \code{"all"}).
#' @param param_type Character. Spatial query type: \code{"bBox"} for a bounding-box
#'   query or \code{"countycode"} for a WQP county-code query.
#' @param param_value Character. The spatial parameter value — a comma-separated
#'   bounding box string (\code{"west,south,east,north"}) or a WQP county code
#'   (\code{"US:12:057"}).
#' @param char_names Character vector. WQP \code{characteristicName} values to
#'   request; passed as repeated query parameters via \code{.multi = "explode"}.
#' @param start_date Character. Query start date in \code{"MM-DD-YYYY"} format.
#' @param end_date Character. Query end date in \code{"MM-DD-YYYY"} format.
#' @param cache_dir Character. Directory where the zip file is saved.
#' @param verbose Logical. Print download progress and file size. Default \code{TRUE}.
#'
#' @return Character path to the zip file on success, or \code{NULL} if the
#'   download fails after retries.

wqp_download_chunk <- function(
  endpoint,
  label,
  param_type,
  param_value,
  char_names,
  start_date,
  end_date,
  cache_dir,
  data_profile = "resultPhysChem",
  verbose = TRUE
) {
  zip_path <- file.path(
    cache_dir,
    sprintf("wqp_%s_%s.zip", tolower(endpoint), label)
  )

  if (file.exists(zip_path)) {
    if (verbose)
      cat(sprintf("  [%s] Using cached %s file\n", label, endpoint))
    return(zip_path)
  }

  if (verbose) cat(sprintf("  [%s] Downloading %s...\n", label, endpoint))

  req <- httr2::request(
    sprintf("https://www.waterqualitydata.us/data/%s/search", endpoint)
  ) |>
    httr2::req_url_query(
      characteristicName = char_names,
      startDateLo        = start_date,
      startDateHi        = end_date,
      mimeType           = "csv",
      zip                = "yes",
      .multi             = "explode"
    ) |>
    httr2::req_timeout(600) |>
    httr2::req_retry(max_tries = 3, backoff = ~ 30)

  # Skipping sort improves throughput on large Result downloads
  if (endpoint == "Result") {
    req <- req |> httr2::req_url_query(sorted = "no")
  }

  if (!is.null(data_profile)) {
    req <- req |> httr2::req_url_query(dataProfile = data_profile)
  }

  # Add spatial param last so I() is not re-encoded by subsequent req_url_query calls
  if (param_type == "bBox") {
    req <- req |> httr2::req_url_query(bBox = param_value)
  } else {
    req <- req |> httr2::req_url_query(countycode = I(param_value))
  }

  if (verbose) message("  [DEBUG] URL: ", req$url)

  tryCatch(
    {
      resp <- httr2::req_perform(req)
      writeBin(httr2::resp_body_raw(resp), zip_path)
      if (verbose) {
        size_mb <- file.size(zip_path) / 1e6
        cat(sprintf("  [%s] %s saved (%.1f MB)\n", label, endpoint, size_mb))
      }
      zip_path
    },
    error = function(e) {
      body <- tryCatch(
        httr2::resp_body_string(e$response),
        error = function(e2) "(no body)"
      )
      message(sprintf("  [%s] %s download failed: %s\n  WQP response: %s",
        label, endpoint, e$message, body))
      if (file.exists(zip_path)) file.remove(zip_path)
      NULL
    }
  )
}

#' Read selected columns from the first CSV inside a WQP zip archive
#'
#' Extracts the first file from a zip archive returned by the Water Quality
#' Portal into \code{cache_dir}, reads the requested columns with
#' \code{readr::read_csv()}, normalises WQP's slash-delimited column names
#' (e.g. \code{ResultMeasure/MeasureUnitCode} becomes
#' \code{ResultMeasure_MeasureUnitCode}), then deletes the extracted CSV before
#' returning.
#'
#' @param zip_path Character. Path to the WQP zip file.
#' @param col_select Character vector. Column names to retain after name
#'   normalisation. Unrecognised names are silently ignored via
#'   \code{dplyr::any_of()}.
#' @param cache_dir Character. Directory used as the extraction destination;
#'   the extracted CSV is removed before the function returns.
#'
#' @return A \code{tibble} containing only the columns in \code{col_select}
#'   that were present in the CSV.

wqp_read_zip_csv <- function(zip_path, col_select, cache_dir) {
  csv_name <- utils::unzip(zip_path, list = TRUE)$Name[1]
  utils::unzip(zip_path, files = csv_name, exdir = cache_dir, overwrite = TRUE)
  csv_path <- file.path(cache_dir, csv_name)
  on.exit(if (file.exists(csv_path)) file.remove(csv_path))

  df <- readr::read_csv(csv_path, show_col_types = FALSE, name_repair = "minimal",
    col_types = readr::cols(ResultMeasureValue = readr::col_character()))
  names(df) <- gsub("/", "_", names(df), fixed = TRUE)
  dplyr::select(df, dplyr::any_of(col_select))
}

#' Fetch salinity observations from the USEPA Water Quality Portal
#'
#' Queries the WQP Result and Station endpoints for salinity measurements within
#' the bounding box of \code{counties} (default) or county by county using WQP
#' county codes (\code{by_county = TRUE}). Downloads are cached to \code{cache_dir}
#' so re-running skips completed chunks. Returns long-term mean salinity per
#' station as an \code{sf} point layer.
#'
#' Units "psu", "PSU", "ppt", and "ppth" are treated as equivalent (~1:1 at
#' estuarine salinities) and pooled before computing station means.
#'
#' @param counties An \code{sf} polygon of the study area. Must have a \code{county}
#'   column when \code{by_county = TRUE}.
#' @param start_date Character. Query start date in \code{"MM-DD-YYYY"} format.
#'   Defaults to five years before today.
#' @param end_date Character. Query end date in \code{"MM-DD-YYYY"} format.
#'   Defaults to today.
#' @param char_names Character vector. WQP characteristic names to request.
#'   Default \code{c("Salinity", "Salinity, water")}.
#' @param valid_units Character vector. Unit codes to retain (case-insensitive).
#'   Default \code{c("psu", "PSU", "ppt", "ppth")}.
#' @param by_county Logical. If \code{FALSE} (default), issues a single bounding-box
#'   query. If \code{TRUE}, loops over counties using WQP \code{countycode} parameters
#'   — use this if the full-extent download is too large or times out.
#' @param cache_dir Character. Directory for cached zip downloads. Created if it
#'   does not exist. Defaults to \code{data-raw/wqp_cache} in the project root.
#' @param crs Integer EPSG code for the output CRS. Default \code{3087}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return An \code{sf} point layer in \code{crs} with columns
#'   \code{MonitoringLocationIdentifier}, \code{salinity_mean} (mean PSU/PPT
#'   across the query period), and \code{salinity_n} (observation count).

fetch_wqp_salinity <- function(
  counties,
  start_date  = format(Sys.Date() - round(5 * 365.25), "%m-%d-%Y"),
  end_date    = format(Sys.Date(), "%m-%d-%Y"),
  char_names  = "Salinity",
  valid_units = c("psu", "PSU", "ppt", "ppth"),
  by_county   = FALSE,
  cache_dir   = here::here("data-raw", "wqp_cache"),
  crs         = 3087L,
  verbose     = TRUE
) {

  # WQP county codes for the seven TBCMP counties (Florida state FIPS = 12)
  county_fips <- c(
    Citrus       = "US:12:017",
    Hernando     = "US:12:053",
    Hillsborough = "US:12:057",
    Manatee      = "US:12:081",
    Pasco        = "US:12:101",
    Pinellas     = "US:12:103",
    Sarasota     = "US:12:115"
  )

  if (!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)

  counties_4326 <- sf::st_transform(counties, 4326)

  # --- build download chunks -----------------------------------------------
  # Each chunk is a named character vector: names = labels, values = spatial
  # parameter values. param_type selects the WQP query parameter name.

  if (by_county) {
    if (!"county" %in% names(counties))
      stop("'counties' must have a 'county' column when by_county = TRUE")
    bad <- setdiff(counties$county, names(county_fips))
    if (length(bad))
      stop("No WQP county-code mapping for: ", paste(bad, collapse = ", "))
    chunks     <- county_fips[counties$county]
    param_type <- "countycode"
  } else {
    bb <- sf::st_bbox(counties_4326)
    chunks <- c(all = paste(
      round(c(bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"]), 5),
      collapse = ","
    ))
    param_type <- "bBox"
  }

  # --- download loop -------------------------------------------------------

  result_list <- vector("list", length(chunks))

  for (i in seq_along(chunks)) {
    label <- names(chunks)[i]
    if (verbose) cat(sprintf("\nChunk %d / %d: %s\n", i, length(chunks), label))

    res_zip <- wqp_download_chunk(
      "Result", label, param_type, chunks[[i]],
      char_names, start_date, end_date, cache_dir,
      data_profile = "resultPhysChem", verbose
    )

    if (!is.null(res_zip)) {
      result_list[[i]] <- wqp_read_zip_csv(res_zip, c(
        "MonitoringLocationIdentifier",
        "LongitudeMeasure",
        "LatitudeMeasure",
        "ActivityStartDate",
        "CharacteristicName",
        "ResultMeasureValue",
        "ResultMeasure_MeasureUnitCode"
      ), cache_dir)
    }
  }

  # --- combine and clean ---------------------------------------------------

  results <- dplyr::bind_rows(result_list) |>
    dplyr::mutate(
      salinity = suppressWarnings(as.numeric(ResultMeasureValue)),
      unit     = tolower(trimws(ResultMeasure_MeasureUnitCode)),
      lon      = suppressWarnings(as.numeric(LongitudeMeasure)),
      lat      = suppressWarnings(as.numeric(LatitudeMeasure))
    ) |>
    dplyr::filter(
      unit %in% tolower(valid_units),
      !is.na(salinity),
      dplyr::between(salinity, 0, 45),  # physical plausibility check
      !is.na(lon), !is.na(lat)
    )

  if (verbose) {
    cat(sprintf("\nObservations after unit/range filter: %d\n", nrow(results)))
    cat("Unit breakdown:\n")
    print(dplyr::count(results, unit, sort = TRUE))
  }

  # --- mean per station, return sf -----------------------------------------

  sal_sf <- results |>
    dplyr::group_by(MonitoringLocationIdentifier) |>
    dplyr::summarise(
      salinity_mean = mean(salinity, na.rm = TRUE),
      salinity_n    = dplyr::n(),
      lon           = dplyr::first(lon),
      lat           = dplyr::first(lat),
      .groups       = "drop"
    ) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
    sf::st_transform(crs)

  if (verbose)
    cat(sprintf("Final spatial layer: %d stations\n", nrow(sal_sf)))

  sal_sf
}
