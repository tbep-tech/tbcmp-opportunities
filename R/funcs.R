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
#'   it does not exist. Defaults to \code{data-raw/wqp_cache} in the project
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
  cache_dir = here::here("data-raw", "wqp_cache"),
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
#'   not exist. Defaults to \code{data-raw/fnai} in the project root.
#' @param crs Integer EPSG code for the output CRS. Default \code{3087L}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return An \code{sf} object clipped to \code{cnt} in \code{crs}.

fetch_fnai <- function(
  url,
  cnt,
  cache_dir = here::here("data-raw", "fnai"),
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
#'   \code{here::here("data-raw", "fnai")}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return An \code{sf} polygon object in \code{crs}.

fetch_clip <- function(
  cnt,
  url = 'https://www.fnai.org/shapefiles/CLIP_v4_02.zip',
  max_priority = 3L,
  crs = 3087L,
  cache_dir = here::here('data-raw', 'fnai'),
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
