library(sf)
library(terra)
library(httr2)
library(dplyr)
library(tidyr)
library(fields)
library(elevatr)
library(here)
library(leaflet)

# ============================================================
# 1. DEFINE STUDY AREA
# ============================================================

load(file = here('data/tbcmp_cnt.RData'))

# Get county boundaries — keep 4326 for VDatum, make 3087 version for analysis
tbcmp_cnt_4326 <- tbcmp_cnt |>
  st_transform(4326)

bbox_4326 <- st_bbox(tbcmp_cnt_4326)

tbcmp_cnt_3087 <- tbcmp_cnt

# ============================================================
# 2. FETCH DEM VIA elevatr
# ============================================================
# elevatr::get_elev_raster() fetches AWS Terrain Tiles (NAVD88, meters)
# z = zoom level controls resolution:
#   z = 10 ~ 76m  |  z = 11 ~ 38m  |  z = 12 ~ 19m  |  z = 13 ~ 10m  |  z = 14 ~ 5m
# For a 5-county extent, z = 12 or 13 is a good balance of resolution vs. download size

cat("Downloading DEM via elevatr...\n")

dem_raw <- get_elev_raster(
  locations = tbcmp_cnt_4326,
  z = 12,
  src = "aws",
  clip = "locations" # clips to your polygon automatically
)

# Convert from raster to SpatRaster
dem_raw <- rast(dem_raw)

cat(sprintf(
  "DEM resolution: %.1f x %.1f m (approx)\n",
  res(project(dem_raw, "EPSG:3087"))[1],
  res(project(dem_raw, "EPSG:3087"))[2]
))

cat(sprintf(
  "DEM value range: %.2f to %.2f m\n",
  global(dem_raw, "min", na.rm = TRUE)[[1]],
  global(dem_raw, "max", na.rm = TRUE)[[1]]
))

# Project to 3087
dem_3087 <- project(dem_raw, "EPSG:3087", method = "bilinear")

# Clip/mask — use the 3087 version of your counties layer
dem_3087 <- crop(dem_3087, vect(tbcmp_cnt_3087)) |>
  mask(vect(tbcmp_cnt_3087))

# Confirm it worked
cat(sprintf(
  "Masked DEM range: %.2f to %.2f m\n",
  global(dem_3087, "min", na.rm = TRUE)[[1]],
  global(dem_3087, "max", na.rm = TRUE)[[1]]
))

plot(
  dem_3087,
  main = "DEM (m NAVD88) — Multi-county",
  col = hcl.colors(50, "Terrain"),
  range = c(-2, 10)
)
plot(st_geometry(tbcmp_cnt_3087), add = TRUE, border = "gray30", lwd = 0.5)

# Upper bound: 5 ft = 1.524 m NAVD88
upper_bound_m <- 5 * 0.3048 # 1.524 m

cat(sprintf(
  "\nDEM elevation range: %.2f to %.2f m\n",
  global(dem_3087, "min", na.rm = TRUE)[[1]],
  global(dem_3087, "max", na.rm = TRUE)[[1]]
))

# ============================================================
# 3. BUILD VDatum SAMPLE GRID (in 4326 for API)
# ============================================================

grid_spacing <- 0.05 # degrees (~5 km)

grid_pts <- expand.grid(
  lon = seq(bbox_4326["xmin"], bbox_4326["xmax"], by = grid_spacing),
  lat = seq(bbox_4326["ymin"], bbox_4326["ymax"], by = grid_spacing)
) |>
  as_tibble()

# Clip grid to study area only
grid_sf <- st_as_sf(grid_pts, coords = c("lon", "lat"), crs = 4326)
in_area <- st_intersects(grid_sf, st_union(tbcmp_cnt_4326), sparse = FALSE)[, 1]
grid_pts <- grid_pts[in_area, ]

cat(sprintf("\nQuerying VDatum API for %d points...\n", nrow(grid_pts)))

# ============================================================
# 4. VDatum API QUERY FUNCTION
# ============================================================
# Convert 0 m MLLW -> NAVD88 at each point.
# The returned t_v_value is the NAVD88 elevation of MLLW (a negative number).

query_vdatum <- function(lon, lat, unit = "m") {
  resp <- tryCatch(
    {
      request("https://vdatum.noaa.gov/vdatumweb/api/convert") |>
        req_url_query(
          lon = round(lon, 6),
          lat = round(lat, 6),
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
        req_timeout(15) |>
        req_retry(max_tries = 3, backoff = ~2) |>
        req_perform()
    },
    error = function(e) {
      message("Request error at (", lon, ", ", lat, "): ", e$message)
      NULL
    }
  )

  if (is.null(resp)) {
    return(NA_real_)
  }

  result <- tryCatch(resp_body_json(resp), error = function(e) NULL)
  if (is.null(result)) {
    return(NA_real_)
  }

  # Value is in t_z (not t_v_value) and returned as a string
  val <- suppressWarnings(as.numeric(result$t_z))

  # -999999 is NOAA's sentinel for out-of-coverage
  if (!is.na(val) && val != -999999) {
    return(val)
  }

  NA_real_
}

# ============================================================
# 5. BATCH QUERY WITH RATE LIMITING
# ============================================================

mllw_offsets <- vector("numeric", nrow(grid_pts))

for (i in seq_len(nrow(grid_pts))) {
  mllw_offsets[i] <- query_vdatum(grid_pts$lon[i], grid_pts$lat[i])

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

grid_pts$mllw_navd88_m <- mllw_offsets

n_failed <- sum(is.na(grid_pts$mllw_navd88_m))
cat(sprintf("\n%d points returned NA (offshore/outside coverage)\n", n_failed))

# ============================================================
# 6. INTERPOLATE OFFSET SURFACE
# ============================================================

pts_clean <- drop_na(grid_pts, mllw_navd88_m)

cat(sprintf("Interpolating from %d valid points...\n", nrow(pts_clean)))
cat(sprintf(
  "MLLW offset range: %.3f to %.3f m NAVD88\n",
  min(pts_clean$mllw_navd88_m),
  max(pts_clean$mllw_navd88_m)
))

# Thin-plate spline — appropriate for smooth tidal datum variation
tps_model <- Tps(
  x = as.matrix(pts_clean[, c("lon", "lat")]),
  Y = pts_clean$mllw_navd88_m
)

# Predict onto a fine grid in 4326 (project to 3087 after)
interp_res <- 0.005 # ~500 m at this latitude
lon_seq <- seq(bbox_4326["xmin"], bbox_4326["xmax"], by = interp_res)
lat_seq <- seq(bbox_4326["ymin"], bbox_4326["ymax"], by = interp_res)
pred_grid <- expand.grid(lon = lon_seq, lat = lat_seq)
pred_grid$mllw_offset_m <- predict(tps_model, as.matrix(pred_grid))

# Build raster in 4326, then project to 3087
offset_rast_4326 <- rast(pred_grid, type = "xyz", crs = "EPSG:4326")
offset_rast_3087 <- project(offset_rast_4326, "EPSG:3087", method = "bilinear")

# Mask to study area
offset_rast_3087 <- mask(offset_rast_3087, vect(tbcmp_cnt_3087))

# ============================================================
# 7. ALIGN OFFSET SURFACE TO DEM
# ============================================================
# Resample offset raster to exactly match DEM grid before comparison

mllw_surface <- resample(offset_rast_3087, dem_3087, method = "bilinear")

# Sanity check — grids must align
stopifnot(
  compareGeom(dem_3087, mllw_surface, stopOnError = FALSE)
)

# ============================================================
# 8. CREATE COASTAL STRATUM MASK
# ============================================================

cat("\nCreating coastal stratum layer...\n")

# Stratum: MLLW (spatially varying) to 5 ft (1.524 m) NAVD88
coastal_stratum <- dem_3087 >= mllw_surface & dem_3087 <= upper_bound_m

# Set cells outside stratum to NA
coastal_stratum[!coastal_stratum] <- NA

# ============================================================
# 9. EXPORT
# ============================================================

# Vector polygon version
cat("Converting to vector...\n")
stratum_poly <- as.polygons(coastal_stratum, dissolve = TRUE) |>
  st_as_sf() |>
  st_make_valid() |>
  mutate(
    stratum = "coastal",
    lower_datum = "MLLW (spatially varying, NAVD88)",
    upper_elev_ft = 5,
    upper_elev_m = upper_bound_m,
    crs = "EPSG:3087"
  )

st_write(stratum_poly, "coastal_stratum_3087.gpkg", delete_dsn = TRUE)

# make leaflet plot for stratum_poly and original Tampa Bay coastal stratum for comparison

# original TB
load(file = '../hmpu-workflow/data/coastal.RData')
coastal_4326 <- st_transform(coastal, 4326)
# drop z
coastal_4326 <- st_zm(coastal_4326)
stratum_poly_4326 <- st_transform(stratum_poly, 4326)


leaflet() |>
  addTiles() |>
  addPolygons(
    data = coastal_4326,
    color = "blue",
    fillOpacity = 0.5,
    group = "Original TB Stratum"
  ) |>
  addPolygons(
    data = stratum_poly_4326,
    color = "red",
    fillOpacity = 0.5,
    group = "Interpolated Stratum"
  ) |>
  addLayersControl(
    overlayGroups = c("Original TB Stratum", "Interpolated Stratum"),
    options = layersControlOptions(collapsed = FALSE)
  )
