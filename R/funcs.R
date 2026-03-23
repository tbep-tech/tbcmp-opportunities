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
