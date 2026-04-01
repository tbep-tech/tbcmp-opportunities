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

  if (is.null(resp)) {
    return(NA_real_)
  }

  result <- tryCatch(httr2::resp_body_json(resp), error = function(e) NULL)
  if (is.null(result)) {
    return(NA_real_)
  }

  val <- suppressWarnings(as.numeric(result$t_z))

  if (!is.na(val) && val != -999999) {
    return(val)
  }

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
  in_area <- sf::st_intersects(
    grid_sf,
    sf::st_union(counties_4326),
    sparse = FALSE
  )[, 1]
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
      cat(sprintf(
        "  %d / %d complete (%.0f%%)\n",
        i,
        nrow(grid_pts),
        100 * i / nrow(grid_pts)
      ))
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
  offset_rast_3087 <- terra::project(
    offset_rast_4326,
    "EPSG:3087",
    method = "bilinear"
  )
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

make_coastal_stratum <- function(
  dem,
  mllw_surface,
  upper_bound_m = 5 * 0.3048
) {
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
    'Sellers',
    'Zephyr',
    'Lacoochee',
    'Okeelanta',
    'Terra Ceia',
    'Weekiwachee',
    'Samsula',
    'Homosassa',
    'Tobiska',
    'Bessie',
    'Delray',
    'Floridana',
    'Kesson',
    'Wulfert',
    'Manatee',
    'Okeechobee',
    'Palmetto',
    'Broward',
    'Matlacha',
    'Haplaquents',
    'Hydraquents',
    'Holopaw',
    'Nittaw',
    'Copeland',
    'Waccasassa'
  )

  series_mesic <- c(
    'Wauchula',
    'Pomona',
    'Pineda',
    'Felda',
    'Myakka',
    'Ora',
    'Vero',
    'Immokalee',
    'Smyrna',
    'Basinger',
    'Anclote',
    'Pompano',
    'EauGallie',
    'Eau Gallie',
    'Chobee',
    'Paisley',
    'Arredondo',
    'Cassia',
    'Blichton',
    'Flemington',
    'Placid',
    'Aripeka',
    'Wabasso',
    'Ona',
    'Pinellas',
    'St. Johns',
    'Winder',
    'Bradenton',
    'Canova',
    'Malabar',
    'Seffner',
    'Parkwood',
    'Braden',
    'Eaton',
    'Farmton',
    'Kanapaha',
    'Lynne',
    'Nobleton',
    'Oldsmar',
    'Waveland',
    'Brynwood',
    'Jumper',
    'Lutterloh',
    'Mabel',
    'Masaryk',
    'Pople',
    'Punta',
    'Tarrytown',
    'Wekiva'
  )

  series_xeric <- c(
    'Tavares',
    'Sparr',
    'Adamsville',
    'Astatula',
    'Chandler',
    'Electra',
    'Paola',
    'Narcoossee',
    'Lake',
    'Pomello',
    'Kendrick',
    'Lochloosa',
    'Newnan',
    'Gainesville',
    'Micanopy',
    'Millhopper',
    'Orlando',
    'Zofo',
    'Zolfo',
    'Candler',
    'Nobleston',
    'Beaches',
    'Pits',
    'Urban land',
    'Quartzipsamments',
    'Arents',
    'Dumps',
    'Canaveral',
    'Orsino',
    'Palm Beach',
    'St. Augustine',
    'Apopka',
    'Archbold',
    'Duette',
    'Fort Meade',
    'Gypsum land',
    'Jonesville',
    'Neilhurst',
    'Udorthents',
    'Slickens',
    'Citronelle',
    'Florahome',
    'Ft. Green',
    'Jonathan',
    'Redlevel',
    'Shadeville',
    'St. Lucie',
    'Sumterville',
    'Williston'
  )

  pat_hydric <- paste0('\\b(', paste(series_hydric, collapse = '|'), ')\\b')
  pat_mesic <- paste0('\\b(', paste(series_mesic, collapse = '|'), ')\\b')
  pat_xeric <- paste0('\\b(', paste(series_xeric, collapse = '|'), ')\\b')

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
        stringr::str_detect(
          muname,
          stringr::regex(pat_xeric, ignore_case = TRUE)
        ) ~ 100L,
        stringr::str_detect(
          muname,
          stringr::regex('muck|depressional|mostly wetland', ignore_case = TRUE)
        ) ~ 300L,
        stringr::str_detect(
          muname,
          stringr::regex(pat_hydric, ignore_case = TRUE)
        ) ~ 300L,
        stringr::str_detect(
          muname,
          stringr::regex(pat_mesic, ignore_case = TRUE)
        ) ~ 200L,
        TRUE ~ NA_integer_
      ),
      Descrip = dplyr::case_when(
        gridcode == 100L ~ 'Xeric',
        gridcode == 200L ~ 'Mesic',
        gridcode == 300L ~ 'Hydric',
        TRUE ~ 'Unclassified'
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

  cat(sprintf(
    '  Querying %d tiles (%.2f deg each)...\n',
    nrow(tiles),
    tile_size
  ))

  results <- vector('list', nrow(tiles))

  for (j in seq_len(nrow(tiles))) {
    results[[j]] <- tryCatch(
      soilDB::SDA_spatialQuery(
        geom = tiles[j, ],
        what = 'mupolygon',
        db = 'SSURGO',
        geomIntersection = TRUE
      ),
      error = function(e) {
        message(sprintf(
          '  Tile %d/%d failed: %s',
          j,
          nrow(tiles),
          conditionMessage(e)
        ))
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
#' @param data_profile Character. WQP \code{dataProfile} value appended to Result
#'   queries. Default \code{"resultPhysChem"}. Set to \code{NULL} to omit.
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
    if (verbose) {
      cat(sprintf("  [%s] Using cached %s file\n", label, endpoint))
    }
    return(zip_path)
  }

  if (verbose) {
    cat(sprintf("  [%s] Downloading %s...\n", label, endpoint))
  }

  req <- httr2::request(
    sprintf("https://www.waterqualitydata.us/data/%s/search", endpoint)
  ) |>
    httr2::req_url_query(
      characteristicName = char_names,
      startDateLo = start_date,
      startDateHi = end_date,
      mimeType = "csv",
      zip = "yes",
      .multi = "explode"
    ) |>
    httr2::req_timeout(600) |>
    httr2::req_retry(max_tries = 3, backoff = ~30)

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

  if (verbose) {
    message("  [DEBUG] URL: ", req$url)
  }

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
      message(sprintf(
        "  [%s] %s download failed: %s\n  WQP response: %s",
        label,
        endpoint,
        e$message,
        body
      ))
      if (file.exists(zip_path)) {
        file.remove(zip_path)
      }
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

  df <- readr::read_csv(
    csv_path,
    show_col_types = FALSE,
    name_repair = "minimal",
    col_types = readr::cols(.default = readr::col_character())
  )
  names(df) <- gsub("/", "_", names(df), fixed = TRUE)
  dplyr::select(df, dplyr::any_of(col_select))
}

#' Fetch WQP salinity data and build a three-class salinity polygon layer
#'
#' Queries the USEPA Water Quality Portal for salinity measurements within
#' \code{counties}, computes long-term mean salinity per station, fits a
#' thin-plate spline to interpolate across the study area, reclassifies into
#' three salinity zones, and returns a smoothed polygon layer. Set
#' \code{return_raw = TRUE} to skip interpolation and return the station point
#' layer instead.
#'
#' Salinity classes follow the original TBCMP schema:
#' \itemize{
#'   \item Class 0 — Fresh (< 0.5 psu)
#'   \item Class 1 — Oligohaline/Mesohaline (0.5–18 psu)
#'   \item Class 3 — Polyhaline/Euhaline (>= 18 psu)
#' }
#'
#' The 18 psu cutoff follows Eleuterius and Eleuterius (1979). Units "psu",
#' "PSU", "ppt", and "ppth" are treated as equivalent (~1:1 at estuarine
#' salinities) and pooled before computing station means.
#'
#' @param counties An \code{sf} polygon of the study area. Must have a
#'   \code{county} column when \code{by_county = TRUE}.
#' @param start_date Character. Query start date in \code{"MM-DD-YYYY"} format.
#' @param end_date Character. Query end date in \code{"MM-DD-YYYY"} format.
#' @param char_names Character vector. WQP characteristic names to request.
#'   Default \code{"Salinity"}.
#' @param valid_units Character vector. Unit codes to retain (case-insensitive).
#'   Default \code{c("psu", "PSU", "ppt", "ppth")}.
#' @param by_county Logical. If \code{FALSE} (default), issues a single
#'   bounding-box query. If \code{TRUE}, loops over counties using WQP
#'   \code{countycode} parameters — use this if the full-extent download is
#'   too large or times out.
#' @param cache_dir Character. Directory for cached zip downloads. Created if
#'   it does not exist. Defaults to \code{data-raw/01_inputs/wqp_cache} in the project
#'   root.
#' @param interp_res Numeric. TPS prediction grid spacing in decimal degrees.
#'   Default \code{0.001} (~100 m). Coarser values (e.g. \code{0.005}) reduce
#'   processing time; \code{smooth = TRUE} compensates for the resulting
#'   blockiness.
#' @param smooth Logical. If \code{TRUE} (default), applies topology-aware
#'   polygon simplification via \code{rmapshaper::ms_simplify()} after
#'   vectorization to remove raster staircase artifacts.
#' @param simplify_keep Numeric (0–1). Fraction of vertices to retain during
#'   simplification. Lower values produce smoother polygons. Default
#'   \code{0.05}. Ignored when \code{smooth = FALSE}.
#' @param return_raw Logical. If \code{TRUE}, returns the station \code{sf}
#'   point layer (with \code{salinity_mean} and \code{salinity_n}) instead of
#'   building the polygon layer. Default \code{FALSE}.
#' @param crs Integer EPSG code for the output CRS. Default \code{3087}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return When \code{return_raw = FALSE} (default), an \code{sf} MULTIPOLYGON
#'   in \code{crs} with columns \code{Classes} (0L, 1L, 3L),
#'   \code{Value_Min}, \code{Value_Max}, and \code{Descrip}. When
#'   \code{return_raw = TRUE}, an \code{sf} point layer with columns
#'   \code{MonitoringLocationIdentifier}, \code{salinity_mean}, and
#'   \code{salinity_n}.

build_salinity_layer <- function(
  counties,
  start_date = '01-01-2016',
  end_date = '12-31-2025',
  char_names = "Salinity",
  valid_units = c("psu", "PSU", "ppt", "ppth"),
  by_county = FALSE,
  cache_dir = here::here("data-raw", "01_inputs", "wqp_cache"),
  interp_res = 0.001,
  smooth = TRUE,
  simplify_keep = 0.05,
  return_raw = FALSE,
  crs = 3087L,
  verbose = TRUE
) {
  # --- validate date format (WQP requires MM-DD-YYYY) ----------------------

  date_rx <- "^(0[1-9]|1[0-2])-(0[1-9]|[12][0-9]|3[01])-[0-9]{4}$"
  if (!grepl(date_rx, start_date)) {
    stop("start_date must be in MM-DD-YYYY format (got: '", start_date, "')")
  }
  if (!grepl(date_rx, end_date)) {
    stop("end_date must be in MM-DD-YYYY format (got: '", end_date, "')")
  }

  # --- fetch WQP salinity data ---------------------------------------------

  county_fips <- c(
    Citrus = "US:12:017",
    Hernando = "US:12:053",
    Hillsborough = "US:12:057",
    Manatee = "US:12:081",
    Pasco = "US:12:101",
    Pinellas = "US:12:103",
    Sarasota = "US:12:115"
  )

  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  counties_4326 <- sf::st_transform(counties, 4326)

  if (by_county) {
    if (!"county" %in% names(counties)) {
      stop("'counties' must have a 'county' column when by_county = TRUE")
    }
    bad <- setdiff(counties$county, names(county_fips))
    if (length(bad)) {
      stop("No WQP county-code mapping for: ", paste(bad, collapse = ", "))
    }
    chunks <- county_fips[counties$county]
    param_type <- "countycode"
  } else {
    bb <- sf::st_bbox(counties_4326)
    chunks <- c(
      all = paste(
        round(c(bb["xmin"], bb["ymin"], bb["xmax"], bb["ymax"]), 5),
        collapse = ","
      )
    )
    param_type <- "bBox"
  }

  result_list <- vector("list", length(chunks))

  for (i in seq_along(chunks)) {
    label <- names(chunks)[i]
    if (verbose) {
      cat(sprintf("\nChunk %d / %d: %s\n", i, length(chunks), label))
    }

    res_zip <- wqp_download_chunk(
      "Result",
      label,
      param_type,
      chunks[[i]],
      char_names,
      start_date,
      end_date,
      cache_dir,
      data_profile = "resultPhysChem",
      verbose
    )

    if (!is.null(res_zip)) {
      result_list[[i]] <- wqp_read_zip_csv(
        res_zip,
        c(
          "MonitoringLocationIdentifier",
          "ActivityLocation_LongitudeMeasure",
          "ActivityLocation_LatitudeMeasure",
          "ActivityStartDate",
          "CharacteristicName",
          "ResultMeasureValue",
          "ResultMeasure_MeasureUnitCode"
        ),
        cache_dir
      )
    }
  }

  results <- dplyr::bind_rows(result_list) |>
    dplyr::mutate(
      salinity = suppressWarnings(as.numeric(ResultMeasureValue)),
      unit = tolower(trimws(`ResultMeasure_MeasureUnitCode`)),
      lon = suppressWarnings(as.numeric(`ActivityLocation_LongitudeMeasure`)),
      lat = suppressWarnings(as.numeric(`ActivityLocation_LatitudeMeasure`))
    ) |>
    dplyr::filter(
      unit %in% tolower(valid_units),
      !is.na(salinity),
      dplyr::between(salinity, 0, 45),
      !is.na(lon),
      !is.na(lat)
    )

  if (verbose) {
    cat(sprintf("\nObservations after unit/range filter: %d\n", nrow(results)))
    cat("Unit breakdown:\n")
    print(dplyr::count(results, unit, sort = TRUE))
  }

  sal_pts <- results |>
    dplyr::group_by(MonitoringLocationIdentifier) |>
    dplyr::summarise(
      salinity_mean = mean(salinity, na.rm = TRUE),
      salinity_n = dplyr::n(),
      lon = dplyr::first(lon),
      lat = dplyr::first(lat),
      .groups = "drop"
    ) |>
    sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
    sf::st_transform(crs)

  if (verbose) {
    cat(sprintf("Station layer: %d stations\n", nrow(sal_pts)))
  }

  if (return_raw) {
    return(sal_pts)
  }

  # --- project points to 4326 for TPS fitting (lon/lat) --------------------

  pts_4326 <- sf::st_transform(sal_pts, 4326)
  coords <- sf::st_coordinates(pts_4326)
  salinity <- pts_4326$salinity_mean

  keep <- !is.na(salinity)
  coords <- coords[keep, , drop = FALSE]
  salinity <- salinity[keep]

  # --- prediction grid -----------------------------------------------------

  bbox_4326 <- sf::st_bbox(sf::st_transform(counties, 4326))

  pred_grid <- expand.grid(
    lon = seq(bbox_4326["xmin"], bbox_4326["xmax"], by = interp_res),
    lat = seq(bbox_4326["ymin"], bbox_4326["ymax"], by = interp_res)
  )

  # --- interpolate via TPS -------------------------------------------------

  if (verbose) {
    cat(sprintf("Fitting TPS to %d stations...\n", nrow(coords)))
  }

  tps_model <- fields::Tps(x = coords, Y = salinity)
  pred_grid$salinity <- pmax(
    0,
    as.numeric(predict(tps_model, as.matrix(pred_grid)))
  )

  if (verbose) {
    cat(sprintf(
      "Predicted salinity range: %.3f to %.3f psu\n",
      min(pred_grid$salinity),
      max(pred_grid$salinity)
    ))
  }

  sal_rast_4326 <- terra::rast(pred_grid, type = "xyz", crs = "EPSG:4326")
  sal_rast_3087 <- terra::project(
    sal_rast_4326,
    "EPSG:3087",
    method = "bilinear"
  )

  # --- reclassify into three classes and vectorize -------------------------
  # Class 0 = Fresh (<0.5), Class 1 = 0.5-18, Class 3 = >=18

  if (verbose) {
    cat("Reclassifying and vectorizing...\n")
  }

  rcl_matrix <- matrix(
    c(
      -Inf,
      0.5,
      0,
      0.5,
      18,
      1,
      18,
      Inf,
      3
    ),
    ncol = 3,
    byrow = TRUE
  )

  rcl_rast <- terra::classify(sal_rast_3087, rcl_matrix)
  class_vals <- c(0L, 1L, 3L)
  descrips <- c("Fresh (<0.5)", "0.5-18", ">18")
  poly_list <- vector("list", 3)

  for (k in seq_along(class_vals)) {
    class_mask <- terra::ifel(rcl_rast == class_vals[k], 1L, NA)
    class_sal <- terra::mask(sal_rast_3087, class_mask)
    val_min <- as.numeric(terra::global(class_sal, "min", na.rm = TRUE))
    val_max <- as.numeric(terra::global(class_sal, "max", na.rm = TRUE))

    poly_list[[k]] <- terra::as.polygons(class_mask, dissolve = TRUE) |>
      sf::st_as_sf() |>
      sf::st_make_valid() |>
      sf::st_cast("MULTIPOLYGON") |>
      sf::st_intersection(sf::st_union(counties)) |>
      sf::st_cast("MULTIPOLYGON") |>
      dplyr::summarise(geometry = sf::st_union(geometry)) |>
      dplyr::mutate(
        Classes = class_vals[k],
        Value_Min = val_min,
        Value_Max = val_max,
        Descrip = descrips[k]
      ) |>
      dplyr::select(Classes, Value_Min, Value_Max, Descrip)
  }

  sal_poly <- dplyr::bind_rows(poly_list) |>
    dplyr::arrange(Classes)

  # --- simplify topology-aware (no gaps between zones) ---------------------

  if (smooth) {
    if (verbose) {
      cat("Simplifying polygons...\n")
    }
    sal_poly <- rmapshaper::ms_simplify(
      sal_poly,
      keep = simplify_keep,
      keep_shapes = TRUE
    )
  }

  if (verbose) {
    cat(sprintf("Final layer: %d features\n", nrow(sal_poly)))
  }

  sal_poly
}

#' Download and clip an FNAI zipped GDB layer to the study area
#'
#' Downloads a zipped File Geodatabase from the FNAI website (or any compatible
#' URL), extracts it, reads the first layer, reprojects, fixes geometries, and
#' clips to the study area boundary. The zip is cached in \code{cache_dir} so
#' subsequent calls skip the download.
#'
#' @param url Character. Direct URL to the FNAI zip file.
#' @param cnt An \code{sf} polygon of the study area used to clip the
#'   output.
#' @param cache_dir Character. Directory for the cached zip. Created if it does
#'   not exist. Defaults to \code{data-raw/01_inputs/fnai} in the project root.
#' @param crs Integer EPSG code for the output CRS. Default \code{3087L}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return An \code{sf} object clipped to \code{cnt} in \code{crs}.

fetch_fnai <- function(
  url,
  cnt,
  cache_dir = here::here("data-raw", "01_inputs", "fnai"),
  crs = 3087L,
  verbose = TRUE
) {
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  zip_path <- file.path(cache_dir, basename(url))

  if (file.exists(zip_path)) {
    if (verbose) cat(sprintf("  Using cached: %s\n", basename(url)))
  } else {
    if (verbose) {
      cat(sprintf("  Downloading: %s\n", basename(url)))
    }
    resp <- httr2::request(url) |>
      httr2::req_timeout(600) |>
      httr2::req_retry(max_tries = 3, backoff = ~30) |>
      httr2::req_perform()
    writeBin(httr2::resp_body_raw(resp), zip_path)
    if (verbose) {
      cat(sprintf("  Saved (%.1f MB)\n", file.size(zip_path) / 1e6))
    }
  }

  tmp_dir <- tempfile()
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE))
  utils::unzip(zip_path, exdir = tmp_dir)

  # prefer GDB, fall back to shapefile
  gdb <- list.files(
    tmp_dir,
    pattern = "\\.gdb$",
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = TRUE
  )
  if (length(gdb)) {
    layer <- sf::st_layers(gdb[1])$name[1]
    tmp_gpkg <- tempfile(fileext = ".gpkg")
    sf::gdal_utils(
      util = "vectortranslate",
      source = gdb[1],
      destination = tmp_gpkg,
      options = c(
        "-nlt",
        "CONVERT_TO_LINEAR",
        "-nlt",
        "PROMOTE_TO_MULTI",
        "-f",
        "GPKG",
        "-lco",
        "SPATIAL_INDEX=NO"
      )
    )
    out <- sf::st_read(tmp_gpkg, quiet = !verbose)
  } else {
    shp <- list.files(
      tmp_dir,
      pattern = "\\.shp$",
      full.names = TRUE,
      recursive = TRUE
    )
    if (!length(shp)) {
      stop("No .gdb or .shp found in zip: ", basename(url))
    }
    out <- sf::st_read(shp[1], quiet = !verbose)
  }

  out |>
    sf::st_zm(drop = TRUE) |>
    sf::st_transform(crs) |>
    sf::st_make_valid() |>
    sf::st_buffer(dist = 0) |>
    sf::st_intersection(sf::st_union(cnt)) |>
    sf::st_make_valid()
}

#' Fetch Florida DEP Aquatic Preserves clipped to a county boundary
#'
#' Downloads the statewide Aquatic Preserves GeoJSON from the Florida DEP
#' ArcGIS service, routes it through \code{gdal_utils} vectortranslate to
#' linearize any curve geometries, then clips to \code{cnt}.
#'
#' @param cnt An \code{sf} object defining the clipping boundary.
#' @param url Character. GeoJSON endpoint URL.
#' @param crs Integer EPSG code for the output CRS. Default \code{3087L}.
#'
#' @return An \code{sf} object clipped to \code{cnt} in \code{crs}.

fetch_aqprs <- function(
  cnt,
  url = 'https://geodata.dep.state.fl.us/datasets/81841412d3984e9aac2c00c21e41d32e_0.geojson',
  crs = 3087L
) {
  tmp <- tempfile(fileext = '.gpkg')
  on.exit(unlink(tmp), add = TRUE)

  sf::gdal_utils(
    util = 'vectortranslate',
    source = url,
    destination = tmp,
    options = c(
      '-nlt',
      'PROMOTE_TO_MULTI',
      '-nlt',
      'CONVERT_TO_LINEAR',
      '-f',
      'GPKG',
      '-lco',
      'SPATIAL_INDEX=NO'
    )
  )

  sf::st_read(tmp, quiet = TRUE) |>
    sf::st_transform(crs) |>
    sf::st_make_valid() |>
    sf::st_buffer(dist = 0) |>
    sf::st_intersection(sf::st_union(cnt)) |>
    sf::st_make_valid()
}

#' Download a single county LULC shapefile from SWFWMD
#'
#' Downloads the 2023 Land Use Land Cover shapefile for one TBCMP county from
#' the Southwest Florida Water Management District (SWFWMD) ArcGIS portal,
#' reprojects to EPSG:3087, and returns the result. The temporary zip and
#' extraction directory are deleted before returning.
#'
#' Source: \url{https://swfwmd.maps.arcgis.com}
#'
#' @param county Character. County name — one of \code{"Citrus"},
#'   \code{"Hernando"}, \code{"Hillsborough"}, \code{"Manatee"},
#'   \code{"Pasco"}, \code{"Pinellas"}, or \code{"Sarasota"}.
#'
#' @return An \code{sf} polygon object in EPSG:3087 with a \code{FLUCCSCODE} column.

fetch_lulc <- function(county) {
  lulc_items <- c(
    Citrus = "ef11d576fcb44a44b8985a2bbc38a4f7",
    Hernando = "86997a7802584735a8d8193f4a1af30f",
    Hillsborough = "f95454790b724f7fb4a396c5de3c94d2",
    Manatee = "5b0895fe575e4ab297e50546d99e407a",
    Pasco = "81b1498949cd4ba8b9db52ca9e47653b",
    Pinellas = "f0d97c30ada6487eb91196de088ee2fc",
    Sarasota = "06b95375e3dc48e1b61a9b95a87aba30"
  )

  if (!county %in% names(lulc_items)) {
    stop("county must be one of: ", paste(names(lulc_items), collapse = ", "))
  }

  url <- paste0(
    "https://swfwmd.maps.arcgis.com/sharing/rest/content/items/",
    lulc_items[[county]],
    "/data"
  )

  message("Downloading LULC for ", county, " ...")

  zip_path <- tempfile(fileext = ".zip")
  on.exit(unlink(zip_path), add = TRUE)
  curl::curl_download(url, zip_path)

  tmp_dir <- tempfile()
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)
  unzip(zip_path, exdir = tmp_dir)

  shp <- list.files(
    tmp_dir,
    pattern = "\\.shp$",
    full.names = TRUE,
    recursive = TRUE
  )[1]

  sf::st_read(shp, quiet = TRUE) |>
    sf::st_transform(3087L) |>
    sf::st_make_valid()
}

#' Download and combine SWFWMD seagrass layers
#'
#' Downloads the 2024 seagrass mapping layers for the Suncoast and Springs Coast
#' study areas from the SWFWMD ArcGIS REST services and combines them into a
#' single layer. Both sources are assumed to have no spatial overlap; a
#' \code{source} column is added before combining to preserve provenance.
#'
#' Sources:
#' \itemize{
#'   \item Suncoast: \url{https://data-swfwmd.opendata.arcgis.com/datasets/swfwmd::seagrass-in-2024/about}
#'   \item Springs Coast: \url{https://data-swfwmd.opendata.arcgis.com/datasets/swfwmd::seagrass-in-2024-for-the-springs-coast/about}
#' }
#'
#' @return An \code{sf} polygon object in EPSG:3087 with columns \code{source},
#'   \code{FLUCCSCODE}, and \code{FLUCCSDESC}.

fetch_seagrass <- function() {
  urls <- c(
    suncoast = "https://www45.swfwmd.state.fl.us/arcgis12/rest/services/OpenData/Environmental_Seagrass2018_sql/MapServer/3/query?outFields=*&where=1%3D1&f=geojson",
    springs_coast = "https://www45.swfwmd.state.fl.us/arcgis12/rest/services/OpenData/Env_sg_springscoast/MapServer/4/query?outFields=*&where=1%3D1&f=geojson"
  )

  layers <- lapply(names(urls), function(src) {
    message("Downloading seagrass: ", src, " ...")
    tmp <- tempfile(fileext = ".gpkg")
    on.exit(unlink(tmp), add = TRUE)
    sf::gdal_utils(
      util = "vectortranslate",
      source = urls[[src]],
      destination = tmp,
      options = c(
        "-nlt",
        "PROMOTE_TO_MULTI",
        "-nlt",
        "CONVERT_TO_LINEAR",
        "-f",
        "GPKG",
        "-lco",
        "SPATIAL_INDEX=NO"
      )
    )
    sf::st_read(tmp, quiet = TRUE) |>
      dplyr::mutate(source = src, .before = 1) |>
      dplyr::select(source, FLUCCSCODE, FLUCCSDESC, geom)
  })

  dplyr::bind_rows(layers) |>
    sf::st_transform(3087L) |>
    sf::st_make_valid()
}

#' Clip combined seagrass layer to a single county
#'
#' Subsets the combined seagrass layer returned by \code{fetch_seagrass()} to
#' the boundary of one county.
#'
#' @param seagrass_all An \code{sf} object as returned by \code{fetch_seagrass()}.
#' @param tbcmp_cnt An \code{sf} polygon with one row per county and a
#'   \code{county} column.
#' @param county Character. County name matching a value in \code{tbcmp_cnt$county}.
#'
#' @return An \code{sf} polygon object clipped to the county boundary in EPSG:3087.

clip_seagrass <- function(seagrass_all, tbcmp_cnt, county) {
  cnt_geom <- tbcmp_cnt[tbcmp_cnt$county == county, ]
  sf::st_intersection(seagrass_all, sf::st_union(cnt_geom)) |>
    sf::st_make_valid()
}

#' Fetch FNAI CLIP priorities raster and vector clipped to a county boundary
#'
#' Downloads the CLIP v4 zip from FNAI, extracts the priorities GDB, reads the
#' \code{CLIPprio_v4} raster via the GDAL OpenFileGDB raster driver (requires
#' GDAL >= 3.7), filters to the requested priority levels, masks to \code{cnt},
#' and polygonises the result. The zip is deleted after processing.
#'
#' @param cnt An \code{sf} object defining the clipping/masking boundary.
#' @param url Character. URL to the CLIP v4 zip file.
#' @param max_priority Integer. Retain priorities 1 through this value. Default \code{3L}.
#' @param crs Integer EPSG code for the output CRS. Default \code{3087L}.
#' @param cache_dir Character. Directory for the downloaded zip. Default
#'   \code{here::here("data-raw", "01_inputs", "fnai")}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return An \code{sf} polygon object in \code{crs}.

fetch_clip <- function(
  cnt,
  url = 'https://www.fnai.org/shapefiles/CLIP_v4_02.zip',
  max_priority = 3L,
  crs = 3087L,
  cache_dir = here::here('data-raw', '01_inputs', 'fnai'),
  verbose = TRUE
) {
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  zip_path <- file.path(cache_dir, basename(url))

  if (verbose) {
    cat('Downloading CLIP zip...\n')
  }
  resp <- httr2::request(url) |>
    httr2::req_timeout(1800) |>
    httr2::req_retry(max_tries = 3, backoff = ~60) |>
    httr2::req_perform()
  writeBin(httr2::resp_body_raw(resp), zip_path)
  if (verbose) {
    cat(sprintf('  Saved (%.0f MB)\n', file.size(zip_path) / 1e6))
  }

  tmp_dir <- tempfile()
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  utils::unzip(zip_path, exdir = tmp_dir)
  clip_gdb <- list.files(
    tmp_dir,
    pattern = 'priorities\\.gdb$',
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = TRUE
  )

  if (!length(clip_gdb)) {
    stop('No priorities.gdb found in CLIP zip')
  }
  r <- terra::rast(paste0('OpenFileGDB:"', clip_gdb[1], '":CLIPprio_v4'))
  r <- terra::project(r, paste0('EPSG:', crs))

  # Build reclassification matrix from the raster's actual pixel values.
  # terra::cats(), sf::st_layers(), and gdalinfo all fail to expose the VAT for
  # this OpenFileGDB raster. Instead, derive the mapping from terra::freq():
  # the raster has Mask Flags: PER_DATASET so "No Resource Value" pixels are
  # already NA, leaving only the 5 priority pixel values. By CLIP v4 convention
  # (confirmed in ArcMap), higher pixel value = higher conservation priority,
  # so the maximum pixel value maps to Priority 1 (highest).
  frq <- terra::freq(r)
  frq <- frq[order(frq$value), ] # sort low -> high pixel value
  frq$priority <- rev(seq_len(nrow(frq))) # max pixel value -> priority 1
  rcl <- as.matrix(frq[, c("value", "priority")])
  r <- terra::classify(r, rcl, others = NA)
  r[r > max_priority] <- NA
  r <- terra::mask(r, terra::vect(cnt))

  v <- terra::as.polygons(r) |>
    sf::st_as_sf() |>
    sf::st_make_valid()

  names(v)[1] <- 'priority'
  v <- v[v$priority >= 1, ]
  v$label <- paste('Priority', v$priority)

  # Remove the large zip now that polygons are built
  unlink(zip_path)
  if (verbose) {
    cat('  CLIP zip removed\n')
  }

  v
}

#' Build proposed and existing conservation land layers
#'
#' Combines the individual FNAI/DEP conservation layers into unified existing
#' (\code{exst}) and proposed (\code{prop}) conservation geometries, then
#' corrects overlap so that any area already in \code{exst} is removed from
#' \code{prop}.
#'
#' Existing conservation (\code{exst}) is the union of:
#' \itemize{
#'   \item \code{flma}  – Florida Conservation Lands
#'   \item \code{ffbot} – Future Forever Board of Trustees projects
#'   \item \code{ffa}   – Future Forever acquisitions
#'   \item \code{aqprs} – Florida DEP Aquatic Preserves
#'   \item \code{exstorig} – original TB existing conservation layer downloaded
#'     from the TBEP open data portal
#' }
#'
#' Proposed conservation (\code{prop}) is the union of \code{clip} (CLIP v4
#' priorities 1–3), with any overlap with \code{exst} moved into \code{exst}.
#'
#' @param flma  \code{sf} object. Florida Conservation Lands.
#' @param ffbot \code{sf} object. Future Forever Board of Trustees projects.
#' @param ffa   \code{sf} object. Future Forever acquisitions.
#' @param aqprs \code{sf} object. Florida DEP Aquatic Preserves.
#' @param clip  \code{sf} object. CLIP v4 priority polygons (priorities 1–3).
#' @param exstorig_url Character. URL to the original TB existing conservation
#'   GeoJSON layer.
#' @param crs Integer EPSG code for all outputs. Default \code{3087L}.
#'
#' @return A named list with elements \code{prop} and \code{exst}, each an
#'   \code{sfc_POLYGON} geometry set in \code{crs}.

build_prop_exst <- function(
  flma,
  ffbot,
  ffa,
  aqprs,
  clip,
  exstorig_url = 'https://opendata.arcgis.com/datasets/e977c851f6dc49c48d0729b3cd30cc92_3.geojson',
  crs = 3087L
) {
  # --- existing conservation -----------------------------------------------

  exst <- sf::st_geometry(flma) |>
    sf::st_union() |>
    sf::st_union(sf::st_union(sf::st_geometry(ffbot))) |>
    sf::st_union(sf::st_union(sf::st_geometry(ffa))) |>
    sf::st_union(sf::st_union(sf::st_geometry(aqprs))) |>
    sf::st_buffer(dist = 0)

  # original TB existing conservation layer (routed through gdal_utils to
  # handle any curve geometry in the GeoJSON source)
  exstorig_tmp <- tempfile(fileext = '.gpkg')
  on.exit(unlink(exstorig_tmp), add = TRUE)
  sf::gdal_utils(
    util = 'vectortranslate',
    source = exstorig_url,
    destination = exstorig_tmp,
    options = c(
      '-nlt',
      'PROMOTE_TO_MULTI',
      '-nlt',
      'CONVERT_TO_LINEAR',
      '-f',
      'GPKG',
      '-lco',
      'SPATIAL_INDEX=NO'
    )
  )
  exstorig <- sf::st_read(exstorig_tmp, quiet = TRUE) |>
    sf::st_transform(crs) |>
    sf::st_union() |>
    sf::st_buffer(dist = 0)

  exst <- sf::st_union(exst, exstorig) |>
    sf::st_cast('POLYGON') |>
    sf::st_union() |>
    sf::st_buffer(dist = 0)

  # --- proposed conservation -----------------------------------------------

  prop <- sf::st_geometry(clip) |>
    sf::st_union() |>
    sf::st_buffer(dist = 0)

  # --- correct overlap: anything proposed that overlaps existing -> existing -

  a <- sf::st_set_precision(exst, 1e5)
  b <- sf::st_set_precision(prop, 1e5)

  op1 <- sf::st_difference(a, b) # existing not in proposed
  op2 <- sf::st_difference(b, a) # proposed not in existing
  op3 <- sf::st_intersection(a, b) # overlap -> moves to existing

  prop <- sf::st_cast(op2, 'POLYGON')
  exst <- sf::st_union(op1, op3) |>
    sf::st_geometry() |>
    sf::st_cast('POLYGON')

  list(prop = prop, exst = exst)
}

# ---------------------------------------------------------------------------
# Current layers helpers (ported from hmpu-workflow)
# ---------------------------------------------------------------------------

#' Fix geometries by union, cast to POLYGON, and zero-buffer
#'
#' Unions all features, buffers by zero to resolve topology errors, and casts
#' to individual POLYGON geometries. Returns a bare \code{sfc} so attributes
#' can be re-attached deliberately by the caller.
#'
#' @param dat An \code{sf} or \code{sfc} object.
#' @return An \code{sfc_POLYGON}.

fixgeo <- function(dat) {
  dat |>
    sf::st_union() |>
    sf::st_buffer(dist = 0) |>
    sf::st_geometry() |>
    sf::st_cast('POLYGON') |>
    sf::st_buffer(dist = 0)
}

#' Join FLUCCS codes to a LULC layer and reclassify coastal uplands
#'
#' Joins the FLUCCS lookup table to \code{lulcin} by \code{FLUCCSCODE},
#' identifies native upland polygons that fall within the coastal stratum, and
#' returns them as a separate \code{'Coastal Uplands'} category. The remaining
#' native uplands are trimmed so there is no overlap.
#'
#' @param lulcin  An \code{sf} object with a \code{FLUCCSCODE} column.
#' @param coastal An \code{sf} or \code{sfc} object for the coastal stratum.
#' @param fluccs  A data frame with \code{FLUCCSCODE} and \code{HMPU_TARGETS}
#'   columns.
#' @return An \code{sf} object with a single \code{HMPU_TARGETS} column.

add_coast_up <- function(lulcin, coastal, fluccs) {
  lulc <- lulcin |>
    dplyr::mutate(FLUCCSCODE = as.integer(FLUCCSCODE)) |>
    dplyr::left_join(fluccs, by = 'FLUCCSCODE') |>
    dplyr::select(HMPU_TARGETS)

  uplands <- lulc |>
    dplyr::filter(HMPU_TARGETS == 'Native Uplands') |>
    sf::st_geometry() |>
    sf::st_union() |>
    sf::st_cast('POLYGON')

  coastal_uplands <- uplands |>
    sf::st_intersection(sf::st_union(sf::st_geometry(coastal))) |>
    sf::st_union() |>
    sf::st_cast('POLYGON')
  coastal_uplands <- sf::st_sf(geometry = coastal_uplands) |>
    dplyr::mutate(HMPU_TARGETS = 'Coastal Uplands') |>
    dplyr::select(HMPU_TARGETS) |>
    sf::st_zm()

  if (nrow(coastal_uplands) == 0) {
    return(lulc)
  }

  lulcdiff <- sf::st_difference(
    lulc,
    sf::st_geometry(sf::st_union(sf::st_combine(coastal_uplands)))
  )
  dplyr::bind_rows(lulcdiff, coastal_uplands)
}

#' Build the four current-condition layers for one county
#'
#' Clips the county-wide shared input layers to the specified county boundary,
#' then produces the four spatial outputs used in the TBCMP opportunity analysis:
#'
#' \describe{
#'   \item{\code{nativelyr}}{Existing and proposed conservation lands
#'     containing native (non-restorable) habitat, with \code{HMPU_TARGETS}
#'     and \code{typ} (\code{"Existing"}/\code{"Proposed"}) columns.}
#'   \item{\code{restorelyr}}{Restorable lands within conservation, split into
#'     sub-categories (Native Uplands, Coastal Uplands, Freshwater Wetlands,
#'     Mangrove Forests/Salt Barrens, Salt Marshes) based on soil type, coastal
#'     stratum position, and salinity zone.}
#'   \item{\code{nativersrv}}{Bare \code{sfc}: native habitats in the coastal
#'     stratum that are proposed for conservation or currently unprotected.}
#'   \item{\code{restorersrv}}{Bare \code{sfc}: restorable lands in the coastal
#'     stratum that are proposed for conservation or currently unprotected.}
#' }
#'
#' @param lulc             \code{sf}. County LULC layer with a \code{FLUCCSCODE} column.
#' @param coastal_stratum  \code{sf}. Full study-area coastal stratum; clipped to
#'   \code{county} internally.
#' @param soils            \code{sf}. Full study-area soils layer with \code{gridcode}
#'   column (100 = Xeric, 200/300 = Mesic/Hydric); clipped to \code{county} internally.
#' @param salinity_layer   \code{sf}. Full study-area salinity layer with \code{Descrip}
#'   column; clipped to \code{county} internally.
#' @param prop             \code{sf} or \code{sfc}. Full study-area proposed conservation
#'   lands; clipped to \code{county} internally.
#' @param exst             \code{sf} or \code{sfc}. Full study-area existing conservation
#'   lands; clipped to \code{county} internally.
#' @param fluccs           Data frame. FLUCCS lookup with \code{FLUCCSCODE} and
#'   \code{HMPU_TARGETS} columns.
#' @param tbcmp_cnt        \code{sf} polygon with one row per county and a \code{county}
#'   column, used to derive the clipping boundary.
#' @param county           Character. County name matching a value in \code{tbcmp_cnt$county}.
#'
#' @return A named list with elements \code{nativelyr}, \code{restorelyr},
#'   \code{nativersrv}, and \code{restorersrv}.

build_current_lyrs <- function(
  lulc,
  coastal_stratum,
  soils,
  salinity_layer,
  prop,
  exst,
  fluccs,
  tbcmp_cnt,
  county
) {
  # Clip shared layers to county boundary
  cnt_geom <- sf::st_union(tbcmp_cnt[tbcmp_cnt$county == county, ])
  coastal <- sf::st_intersection(coastal_stratum, cnt_geom)
  soils <- sf::st_intersection(soils, cnt_geom)
  salin <- sf::st_intersection(salinity_layer, cnt_geom)
  prop <- sf::st_intersection(prop, cnt_geom)
  exst <- sf::st_intersection(exst, cnt_geom)

  # Prepare LULC: FLUCCS join, coastal uplands reclassification, drop non-habitat
  lulc_prep <- add_coast_up(lulc, coastal, fluccs) |>
    dplyr::filter(!HMPU_TARGETS %in% c('Developed', 'Open Water'))
  categories <- unique(lulc_prep$HMPU_TARGETS)
  prop_geom <- sf::st_union(prop)
  exst_geom <- sf::st_union(exst)
  coastal_geom <- sf::st_union(coastal)

  # Intersect each HMPU category with proposed and existing conservation
  propall <- NULL
  exstall <- NULL
  for (cat in categories) {
    tmp <- lulc_prep |>
      dplyr::filter(HMPU_TARGETS == cat) |>
      fixgeo()
    propall <- dplyr::bind_rows(
      propall,
      sf::st_sf(geometry = sf::st_intersection(tmp, prop_geom) |> fixgeo()) |>
        dplyr::mutate(HMPU_TARGETS = cat, typ = 'Proposed')
    )
    exstall <- dplyr::bind_rows(
      exstall,
      sf::st_sf(geometry = sf::st_intersection(tmp, exst_geom) |> fixgeo()) |>
        dplyr::mutate(HMPU_TARGETS = cat, typ = 'Existing')
    )
  }

  # Native layer: non-Restorable features in proposed/existing conservation
  nativelyr <- dplyr::bind_rows(propall, exstall) |>
    dplyr::filter(HMPU_TARGETS != 'Restorable')

  # Restorable layer: split by soil type, tidal position, and salinity
  restorable <- dplyr::bind_rows(propall, exstall) |>
    dplyr::filter(HMPU_TARGETS == 'Restorable')
  soilsforest <- soils |> dplyr::filter(gridcode == 100) |> fixgeo()
  soilswetland <- soils |> dplyr::filter(gridcode != 100) |> fixgeo()
  salinlo <- salin |> dplyr::filter(Descrip == '0.5-18') |> fixgeo()

  restorelyr <- NULL
  for (cur_typ in c('Proposed', 'Existing')) {
    tmp <- restorable |>
      dplyr::filter(typ == cur_typ) |>
      fixgeo()

    uplands <- sf::st_intersection(tmp, soilsforest)
    coastal_ups <- sf::st_intersection(uplands, coastal_geom)
    uplands <- sf::st_difference(uplands, sf::st_union(coastal_ups))
    wetlands <- sf::st_intersection(tmp, soilswetland)
    tidal_wetlands <- sf::st_intersection(wetlands, coastal_geom)
    wetlands <- sf::st_difference(wetlands, sf::st_union(tidal_wetlands))
    salt_marshes <- sf::st_intersection(tidal_wetlands, salinlo)
    tidal_wetlands <- sf::st_difference(
      tidal_wetlands,
      sf::st_union(salt_marshes)
    )

    restorelyr <- dplyr::bind_rows(
      restorelyr,
      dplyr::bind_rows(
        sf::st_sf(geometry = fixgeo(uplands)) |>
          dplyr::mutate(HMPU_TARGETS = 'Native Uplands'),
        sf::st_sf(geometry = fixgeo(coastal_ups)) |>
          dplyr::mutate(HMPU_TARGETS = 'Coastal Uplands'),
        sf::st_sf(geometry = fixgeo(wetlands)) |>
          dplyr::mutate(HMPU_TARGETS = 'Freshwater Wetlands'),
        sf::st_sf(geometry = fixgeo(tidal_wetlands)) |>
          dplyr::mutate(HMPU_TARGETS = 'Mangrove Forests/Salt Barrens'),
        sf::st_sf(geometry = fixgeo(salt_marshes)) |>
          dplyr::mutate(HMPU_TARGETS = 'Salt Marshes')
      ) |>
        sf::st_zm(drop = TRUE) |>
        dplyr::mutate(typ = cur_typ)
    )
  }

  # Reservation layers: coastal stratum features proposed or unprotected
  uniexstall <- exstall |> sf::st_union() |> sf::st_make_valid()

  nativersrv <- nativelyr |>
    dplyr::filter(typ == 'Proposed') |>
    sf::st_intersection(coastal_geom) |>
    fixgeo()

  restorersrv <- restorelyr |>
    dplyr::filter(typ == 'Proposed') |>
    sf::st_intersection(coastal_geom) |>
    fixgeo()

  nativeunpro <- lulc_prep |>
    dplyr::filter(HMPU_TARGETS != 'Restorable') |>
    sf::st_intersection(coastal_geom) |>
    fixgeo() |>
    sf::st_difference(uniexstall) |>
    fixgeo()

  restoreunpro <- lulc_prep |>
    dplyr::filter(HMPU_TARGETS == 'Restorable') |>
    sf::st_intersection(coastal_geom) |>
    fixgeo() |>
    sf::st_difference(uniexstall) |>
    fixgeo()

  restorersrv <- c(sf::st_make_valid(restorersrv), restoreunpro) |> fixgeo()
  nativersrv <- c(sf::st_make_valid(nativersrv), nativeunpro) |> fixgeo()

  list(
    nativelyr = nativelyr,
    restorelyr = restorelyr,
    nativersrv = nativersrv,
    restorersrv = restorersrv
  )
}

# ---------------------------------------------------------------------------
# Current extent table helpers (ported from hmpu-workflow)
# ---------------------------------------------------------------------------

#' Estimate LULC area in acres per HMPU target category
#'
#' Filters out subtidal and open-water FLUCCS codes, calls \code{add_coast_up()}
#' to reclassify coastal uplands, then summarises area in acres by
#' \code{HMPU_TARGETS}.
#'
#' @param lulcin An \code{sf} polygon with a \code{FLUCCSCODE} column.
#' @param coastal An \code{sf} or \code{sfc} object for the coastal stratum.
#' @param fluccs A data frame with \code{FLUCCSCODE} and \code{HMPU_TARGETS} columns.
#' @param sumout Logical. If \code{TRUE} (default) return a summary data frame;
#'   if \code{FALSE} return the prepared \code{sf} layer.
#'
#' @return A data frame with columns \code{HMPU_TARGETS} and \code{Acres}, or
#'   the prepared \code{sf} layer when \code{sumout = FALSE}.

lulc_est <- function(lulcin, coastal, fluccs, sumout = TRUE) {
  # FLUCCS codes to remove: subtidal and open-water features tracked separately
  # (bays/estuaries, major water bodies, gulf, tidal flats, oyster bars,
  # submerged sand, patchy/continuous seagrass, attached algae, hardbottom)
  cds <- c(
    5400,
    5700,
    5720,
    6510,
    6540,
    7210,
    9113,
    9116,
    9121,
    9510,
    9511,
    9512,
    9513,
    9514,
    9515
  )

  out <- lulcin %>%
    dplyr::mutate(FLUCCSCODE = as.integer(FLUCCSCODE)) %>%
    dplyr::filter(!FLUCCSCODE %in% cds) %>%
    add_coast_up(coastal, fluccs)

  if (!sumout) {
    return(out)
  }

  out <- out %>%
    dplyr::mutate(
      Acres = sf::st_area(.),
      Acres = units::set_units(Acres, acres),
      Acres = as.numeric(Acres)
    ) %>%
    sf::st_set_geometry(NULL) %>%
    dplyr::group_by(HMPU_TARGETS) %>%
    dplyr::summarise(Acres = sum(Acres), .groups = 'drop') %>%
    dplyr::arrange(HMPU_TARGETS)

  return(out)
}

#' Estimate subtidal area in acres per HMPU target category
#'
#' Joins FLUCCS codes to a seagrass / subtidal layer and summarises area in
#' acres by \code{HMPU_TARGETS}.
#'
#' @param subtin An \code{sf} polygon with a \code{FLUCCSCODE} column.
#' @param fluccs A data frame with \code{FLUCCSCODE} and \code{HMPU_TARGETS} columns.
#' @param sumout Logical. If \code{TRUE} (default) return a summary data frame;
#'   if \code{FALSE} return the prepared \code{sf} layer.
#'
#' @return A data frame with columns \code{HMPU_TARGETS} and \code{Acres}, or
#'   the prepared \code{sf} layer when \code{sumout = FALSE}.

subt_est <- function(subtin, fluccs, sumout = TRUE) {
  out <- subtin %>%
    dplyr::mutate(FLUCCSCODE = as.integer(FLUCCSCODE)) %>%
    dplyr::left_join(fluccs, by = 'FLUCCSCODE') %>%
    dplyr::select(HMPU_TARGETS)

  if (!sumout) {
    return(out)
  }

  out %>%
    dplyr::mutate(
      Acres = sf::st_area(.),
      Acres = units::set_units(Acres, acres),
      Acres = as.numeric(Acres)
    ) %>%
    sf::st_set_geometry(NULL) %>%
    dplyr::group_by(HMPU_TARGETS) %>%
    dplyr::summarise(Acres = sum(Acres), .groups = 'drop') %>%
    dplyr::arrange(HMPU_TARGETS)
}

#' Build a current extent for one spatial unit
#'
#' Clips the coastal stratum to the specified county, then computes current
#' extent, native conservation, and restorable land summaries from the input
#' layers and assembles. Subtidal
#' features (seagrasses, oyster bars) are derived from \code{subt} via
#' \code{subt_est()}. Intertidal and supratidal features are derived from \code{lulc} via \code{lulc_est()}.
#'
#' @param lulc           \code{sf}. County LULC layer with \code{FLUCCSCODE} column.
#' @param subt           \code{sf}. County seagrass / subtidal layer with \code{FLUCCSCODE}.
#' @param coastal_stratum \code{sf}. Full study-area coastal stratum; clipped to
#'   \code{county} internally.
#' @param fluccs         Data frame. FLUCCS lookup with \code{FLUCCSCODE} and
#'   \code{HMPU_TARGETS} columns.
#' @param strata         Data frame. Stratification lookup with \code{Category} and
#'   \code{HMPU_TARGETS} columns (as built in \code{01_inputs.R}).
#' @param nativelyr      \code{sf}. Native habitats in conservation lands (from
#'   \code{build_current_lyrs()}).
#' @param restorelyr     \code{sf}. Restorable lands in conservation (from
#'   \code{build_current_lyrs()}).
#' @param tbcmp_cnt      \code{sf} polygon with one row per county and a \code{county}
#'   column, used to derive the clipping boundary.
#' @param county         Character. County name matching a value in \code{tbcmp_cnt$county}.
#' @param cap            Character. Table caption string.
#'
#' @return A summary of current extent, native conservation, and restorable lands by HMPU
#' target for the county, as a data frame.

curex_fun <- function(
  lulc,
  subt,
  coastal_stratum,
  fluccs,
  strata,
  nativelyr,
  restorelyr,
  tbcmp_cnt,
  county
) {
  # Clip coastal stratum to county boundary
  coastal <- sf::st_intersection(
    coastal_stratum,
    sf::st_union(tbcmp_cnt[tbcmp_cnt$county == county, ])
  )

  # current lulc and subtidal summaries
  lulcsum <- lulc_est(lulc, coastal, fluccs)
  subtsum <- subt_est(subt, fluccs)

  # current summary: join all sources to strata so every target row is present
  cursum <- dplyr::bind_rows(lulcsum, subtsum) %>%
    dplyr::mutate(unis = 'ac', `Current Extent` = Acres) %>%
    dplyr::left_join(strata, ., by = 'HMPU_TARGETS') %>%
    dplyr::filter(!HMPU_TARGETS %in% 'Total Intertidal') %>%
    dplyr::mutate(
      unis = dplyr::if_else(is.na(unis), 'ac', unis),
      `Current Extent` = dplyr::if_else(
        is.na(`Current Extent`),
        0,
        `Current Extent`
      )
    ) %>%
    dplyr::select(Category, HMPU_TARGETS, unis, `Current Extent`) %>%
    dplyr::arrange(Category, HMPU_TARGETS)

  # native habitats in conservation lands summary
  nativesum <- nativelyr %>%
    dplyr::mutate(
      Acres = sf::st_area(.),
      Acres = units::set_units(Acres, acres),
      Acres = as.numeric(Acres),
      typ = paste('native', typ),
      typ = factor(typ, levels = c('native Existing', 'native Proposed'))
    ) %>%
    sf::st_set_geometry(NULL) %>%
    dplyr::group_by(typ, HMPU_TARGETS) %>%
    dplyr::summarise(Acres = sum(Acres), .groups = 'drop') %>%
    dplyr::arrange(typ, HMPU_TARGETS) %>%
    tidyr::spread(typ, Acres, fill = 0, drop = FALSE)

  # restorable lands in conservation summary
  restoresum <- restorelyr %>%
    dplyr::mutate(
      Acres = sf::st_area(.),
      Acres = units::set_units(Acres, acres),
      Acres = as.numeric(Acres),
      typ = paste('restorable', typ),
      typ = factor(
        typ,
        levels = c('restorable Existing', 'restorable Proposed')
      )
    ) %>%
    sf::st_set_geometry(NULL) %>%
    dplyr::group_by(typ, HMPU_TARGETS) %>%
    dplyr::summarise(Acres = sum(Acres, na.rm = TRUE), .groups = 'drop') %>%
    dplyr::arrange(typ, HMPU_TARGETS)

  # split 'Mangrove Forests/Salt Barrens' and 'Freshwater Wetlands' into
  # their two constituent HMPU targets so every target row appears in the table
  duplab1 <- 'Mangrove Forests/Salt Barrens'
  dups1 <- restoresum %>%
    dplyr::filter(HMPU_TARGETS %in% duplab1) %>%
    dplyr::mutate(HMPU_TARGETS = 'Mangrove Forests')

  duplab2 <- 'Freshwater Wetlands'
  dups2 <- restoresum %>%
    dplyr::filter(HMPU_TARGETS %in% duplab2) %>%
    dplyr::mutate(HMPU_TARGETS = 'Non-Forested Freshwater Wetlands')

  restoresum <- restoresum %>%
    dplyr::bind_rows(dups1) %>%
    dplyr::bind_rows(dups2) %>%
    dplyr::mutate(
      HMPU_TARGETS = dplyr::case_when(
        HMPU_TARGETS %in% duplab1 ~ 'Salt Barrens',
        HMPU_TARGETS %in% duplab2 ~ 'Forested Freshwater Wetlands',
        TRUE ~ HMPU_TARGETS
      )
    ) %>%
    tidyr::spread(typ, Acres, fill = 0, drop = FALSE) %>%
    dplyr::mutate(
      `total restorable` = `restorable Existing` + `restorable Proposed`
    )

  out <- cursum %>%
    dplyr::left_join(nativesum, by = 'HMPU_TARGETS') %>%
    dplyr::left_join(restoresum, by = 'HMPU_TARGETS')

  return(out)
}

#' Create current extent table from pre-computed summaries
#'
#' Formats output form \code{curex_fun} into a table using \code{flextable}.
#'
#' @param allsum     Data frame. Current extent summary (output from \code{curex_fun()}).
#' @param county     Character. County name for the caption.
#'
#' @return A \code{flextable} object.

curextab_fun <- function(allsum, county) {
  allsumfrm <- allsum %>%
    tidyr::gather('var', 'val', -Category, -HMPU_TARGETS, -unis) %>%
    dplyr::mutate(
      val = dplyr::case_when(
        !is.na(val) ~ paste(prettyNum(round(val, 0), big.mark = ','), unis),
        TRUE ~ 'N/A'
      ),
      val = dplyr::case_when(
        (HMPU_TARGETS %in% 'Salt Marshes') &
          (var %in%
            c(
              'total restorable',
              'restorable Existing',
              'restorable Proposed'
            )) ~
          paste(val, '(JU)'),
        TRUE ~ val
      ),
      Category = factor(
        Category,
        levels = c('Subtidal', 'Intertidal', 'Supratidal')
      ),
      HMPU_TARGETS = factor(HMPU_TARGETS, levels = levels(strata$HMPU_TARGETS))
    ) %>%
    tidyr::spread(var, val) %>%
    dplyr::select(-unis) %>%
    dplyr::mutate(
      `native Existing` = dplyr::case_when(
        Category == 'Subtidal' ~ `Current Extent`,
        TRUE ~ `native Existing`
      )
    ) %>%
    dplyr::select(
      Category,
      HMPU_TARGETS,
      `Current Extent`,
      `native Existing`,
      `native Proposed`,
      `total restorable`,
      `restorable Existing`,
      `restorable Proposed`
    )

  cap <- paste(
    'Current habitat extent and conservation opportunity -',
    county,
    'County'
  )
  capfrm <- flextable::as_paragraph(
    flextable::as_chunk(
      cap,
      props = flextable::fp_text_default(font.size = 14, bold = TRUE)
    )
  )

  allsumfrm %>%
    flextable::as_grouped_data(groups = 'Category') %>%
    dplyr::mutate(
      HMPU_TARGETS = dplyr::case_when(
        is.na(HMPU_TARGETS) ~ Category,
        TRUE ~ HMPU_TARGETS
      )
    ) %>%
    dplyr::select(-Category) %>%
    flextable::flextable() %>%
    flextable::set_header_labels(
      HMPU_TARGETS = 'Habitat Type',
      `native Existing` = 'Existing Conservation Lands',
      `native Proposed` = 'Proposed Conservation Lands',
      `total restorable` = 'Total Restoration Opportunity',
      `restorable Existing` = 'Existing Conservation Lands Restoration Opportunity',
      `restorable Proposed` = 'Proposed Conservation Lands Restoration Opportunity'
    ) %>%
    flextable::merge_at(i = 1, part = 'body') %>%
    flextable::merge_at(i = 5, part = 'body') %>%
    flextable::merge_at(i = 9, part = 'body') %>%
    flextable::merge_at(i = 6:7, j = 5, part = 'body') %>%
    flextable::merge_at(i = 6:7, j = 6, part = 'body') %>%
    flextable::merge_at(i = 6:7, j = 7, part = 'body') %>%
    flextable::merge_at(i = 11:12, j = 5, part = 'body') %>%
    flextable::merge_at(i = 11:12, j = 6, part = 'body') %>%
    flextable::merge_at(i = 11:12, j = 7, part = 'body') %>%
    flextable::add_header_row(
      colwidths = c(1, 3, 3),
      values = c('', 'Native Habitats', 'Restorable Habitats')
    ) %>%
    flextable::add_footer_lines(values = '') %>%
    flextable::add_footer_lines(
      values = flextable::as_paragraph(
        'N/A - Not Applicable; JU - Potential ',
        flextable::as_i('Juncus'),
        ' Marsh Opportunity'
      )
    ) %>%
    flextable::fontsize(size = 8, part = 'footer') %>%
    flextable::align(align = 'center', part = 'header') %>%
    flextable::align(
      i = c(2:4, 6:8, 10:13),
      j = 2:7,
      align = 'center',
      part = 'body'
    ) %>%
    flextable::bg(i = c(1, 5, 9), bg = '#418DA8', part = 'body') %>%
    flextable::color(i = c(1, 5, 9), color = 'white', part = 'body') %>%
    flextable::bg(i = 1, bg = '#1A596B', part = 'header') %>%
    flextable::bg(i = 2, j = 1, bg = '#1A596B', part = 'header') %>%
    flextable::color(i = 1, color = 'white', part = 'header') %>%
    flextable::color(i = 2, j = 1, color = 'white', part = 'header') %>%
    flextable::border_outer(part = 'body') %>%
    flextable::border_outer(part = 'header') %>%
    flextable::border_inner_h(part = 'body') %>%
    flextable::border_inner_v(part = 'body') %>%
    flextable::border_inner_h(part = 'header') %>%
    flextable::border_inner_v(part = 'header') %>%
    flextable::set_caption(caption = capfrm) %>%
    flextable::font(part = 'all', fontname = 'Roboto')
}

# ---------------------------------------------------------------------------
# Opportunity map helpers
# ---------------------------------------------------------------------------

#' Build the opportunity layer for one county
#'
#' Clips the coastal stratum to the specified county, then combines the four
#' current-condition layers into a single \code{sf} object with a \code{cat}
#' column distinguishing six opportunity categories:
#'
#' \itemize{
#'   \item \code{Reservation Native} — native habitats in the coastal stratum
#'     that are proposed for conservation or currently unprotected
#'   \item \code{Reservation Restorable} — restorable lands in the coastal
#'     stratum that are proposed for conservation or currently unprotected
#'   \item \code{Existing Conservation Native} — native habitats within
#'     existing conservation lands
#'   \item \code{Proposed Conservation Native} — native habitats within
#'     proposed conservation lands outside the coastal stratum
#'   \item \code{Existing Conservation Restorable} — restorable lands within
#'     existing conservation lands
#'   \item \code{Proposed Conservation Restorable} — restorable lands within
#'     proposed conservation lands outside the coastal stratum
#' }
#'
#' Proposed conservation features that overlap the coastal stratum are removed
#' via \code{sf::st_difference()} to avoid double-counting with the reservation
#' layers.
#'
#' @param nativersrv      \code{sfc}. Native habitats in the coastal reservation
#'   space, as returned by \code{build_current_lyrs()}.
#' @param restorersrv     \code{sfc}. Restorable lands in the coastal reservation
#'   space, as returned by \code{build_current_lyrs()}.
#' @param nativelyr       \code{sf}. Native habitats in conservation lands, with
#'   a \code{typ} column (\code{"Existing"} / \code{"Proposed"}).
#' @param restorelyr      \code{sf}. Restorable lands in conservation lands, with
#'   a \code{typ} column.
#' @param coastal_stratum \code{sf}. Full study-area coastal stratum; clipped to
#'   \code{county} internally.
#' @param tbcmp_cnt       \code{sf} polygon with one row per county and a
#'   \code{county} column.
#' @param county          Character. County name matching a value in
#'   \code{tbcmp_cnt$county}.
#'
#' @return An \code{sf} POLYGON object with a single \code{cat} column
#'   identifying the opportunity category.

oppdat_fun <- function(
  nativersrv,
  restorersrv,
  nativelyr,
  restorelyr,
  coastal_stratum,
  tbcmp_cnt,
  county
) {
  # Clip coastal stratum to county and union for differencing
  cnt_geom  <- sf::st_union(tbcmp_cnt[tbcmp_cnt$county == county, ])
  coastal   <- sf::st_intersection(coastal_stratum, cnt_geom)
  unicoastal <- sf::st_union(sf::st_combine(coastal)) |> sf::st_make_valid()

  # Reservation layers (already within the coastal stratum)
  nativersrv_out <- sf::st_sf(geometry = fixgeo(nativersrv)) |>
    dplyr::mutate(cat = 'Reservation Native')

  restorersrv_out <- sf::st_sf(geometry = fixgeo(restorersrv)) |>
    dplyr::mutate(cat = 'Reservation Restorable')

  # Existing conservation layers
  nativelyr_exst <- sf::st_sf(
    geometry = nativelyr |> dplyr::filter(typ == 'Existing') |> fixgeo()
  ) |>
    dplyr::mutate(cat = 'Existing Conservation Native')

  restorelyr_exst <- sf::st_sf(
    geometry = restorelyr |> dplyr::filter(typ == 'Existing') |> fixgeo()
  ) |>
    dplyr::mutate(cat = 'Existing Conservation Restorable')

  # Proposed conservation layers: remove coastal stratum overlap to avoid
  # double-counting with reservation layers
  nativelyr_prop <- nativelyr |>
    dplyr::filter(typ == 'Proposed') |>
    fixgeo() |>
    sf::st_difference(unicoastal) |>
    fixgeo()
  nativelyr_prop <- sf::st_sf(geometry = nativelyr_prop) |>
    dplyr::mutate(cat = 'Proposed Conservation Native') |>
    sf::st_make_valid()

  restorelyr_prop <- restorelyr |>
    dplyr::filter(typ == 'Proposed') |>
    fixgeo() |>
    sf::st_difference(unicoastal) |>
    fixgeo()
  restorelyr_prop <- sf::st_sf(geometry = restorelyr_prop) |>
    dplyr::mutate(cat = 'Proposed Conservation Restorable') |>
    sf::st_make_valid()

  dplyr::bind_rows(
    nativersrv_out,
    restorersrv_out,
    nativelyr_exst,
    nativelyr_prop,
    restorelyr_exst,
    restorelyr_prop
  )
}

#' Interactive leaflet map of opportunity layers
#'
#' Maps the output of \code{oppdat_fun()} using \code{leaflet}, colouring
#' polygons by the six opportunity categories. The colour palette matches
#' the ggplot2 static maps from the original hmpu-workflow:
#'
#' \itemize{
#'   \item Existing Conservation Native — \code{yellowgreen}
#'   \item Existing Conservation Restorable — \code{green4}
#'   \item Proposed Conservation Native — \code{dodgerblue1}
#'   \item Proposed Conservation Restorable — \code{dodgerblue4}
#'   \item Reservation Native — \code{violetred1}
#'   \item Reservation Restorable — \code{violetred3}
#' }
#'
#' @param oppdat An \code{sf} object as returned by \code{oppdat_fun()}, with
#'   a \code{cat} column identifying the opportunity category.
#'
#' @return A \code{leaflet} map object.

oppmap_leaflet <- function(oppdat) {
  cols <- c(
    'Existing Conservation Native'     = 'yellowgreen',
    'Existing Conservation Restorable' = 'green4',
    'Proposed Conservation Native'     = 'dodgerblue1',
    'Proposed Conservation Restorable' = 'dodgerblue4',
    'Reservation Native'               = 'violetred1',
    'Reservation Restorable'           = 'violetred3'
  )

  pal <- leaflet::colorFactor(
    palette = unname(cols),
    levels  = names(cols),
    ordered = TRUE
  )

  oppdat_4326 <- sf::st_transform(oppdat, 4326)

  leaflet::leaflet() |>
    leaflet::addProviderTiles(leaflet::providers$CartoDB.Positron) |>
    leaflet::addPolygons(
      data        = oppdat_4326,
      fillColor   = ~pal(cat),
      fillOpacity = 0.7,
      color       = NA,
      weight      = 0,
      label       = ~cat
    ) |>
    leaflet::addLegend(
      pal      = pal,
      values   = names(cols),
      title    = 'Opportunity',
      position = 'bottomright'
    )
}
