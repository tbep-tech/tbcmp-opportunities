# setup ------------------------------------------------------------------

library(tidyverse)
library(sf)
library(curl)
library(here)
library(terra)
library(httr2)
library(fields)
library(leaflet)
library(soilDB)

source(here('R/funcs.R'))

prj <- 3087

# FLUCCS lookup ----------------------------------------------------------

fluccs <- read.csv(
  here('data-raw', 'FLUCCShabsclass.csv'),
  stringsAsFactors = F
)

save(fluccs, file = here('data', 'fluccs.RData'), compress = 'xz')

# stratification lookup table ---------------------------------------------

strata <- data.frame(
  Category = c(
    "Subtidal",
    "Subtidal",
    "Subtidal",
    "Subtidal",
    "Subtidal",
    "Intertidal",
    "Intertidal",
    "Intertidal",
    "Intertidal",
    "Intertidal",
    "Intertidal",
    "Supratidal",
    "Supratidal",
    "Supratidal",
    "Supratidal"
  ),
  HMPU_TARGETS = c(
    "Hard Bottom",
    "Artificial Reefs",
    "Tidal Flats",
    "Seagrasses",
    "Oyster Bars",
    "Total Intertidal",
    "Mangrove Forests",
    "Salt Barrens",
    "Salt Marshes",
    "Living Shorelines",
    "Tidal Tributaries",
    "Coastal Uplands",
    "Non-Forested Freshwater Wetlands",
    "Forested Freshwater Wetlands",
    "Native Uplands"
  ),
  stringsAsFactors = F
) |>
  mutate(
    Category = factor(
      Category,
      levels = c("Subtidal", "Intertidal", "Supratidal")
    ),
    HMPU_TARGETS = factor(HMPU_TARGETS, levels = HMPU_TARGETS)
  )

save(strata, file = here('data', 'strata.RData'), compress = 'xz')

# tbcmp counties ---------------------------------------------------------

# # TBCMP project area bounded by 7 county census areas including associated coastal areas tied to local/state jurisdictions
# curl_download(url = "https://www2.census.gov/geo/tiger/TIGER2025/COUSUB/tl_2025_12_cousub.zip",
#              destfile = "./data-raw/tl_2025_12_cousub.zip")
# unzip("./data-raw/tl_2025_12_cousub.zip", exdir = "./data-raw")

fl_counties <- st_read("./data-raw/tl_2025_12_cousub.shp")

tbcmp_cnt <- fl_counties |>
  st_transform(prj) |>
  filter(COUNTYFP %in% c("017", "053", "057", "081", "101", "103", "115")) |>
  mutate(
    county = case_when(
      COUNTYFP == "017" ~ "Citrus",
      COUNTYFP == "053" ~ "Hernando",
      COUNTYFP == "101" ~ "Pasco",
      COUNTYFP == "103" ~ "Pinellas",
      COUNTYFP == "057" ~ "Hillsborough",
      COUNTYFP == "081" ~ "Manatee",
      COUNTYFP == "115" ~ "Sarasota",
      TRUE ~ "MISSING"
    )
  ) |>
  group_by(county) |>
  summarise() |>
  ungroup()

save(tbcmp_cnt, file = here('data', 'tbcmp_cnt.RData'), compress = 'xz')

# coastal stratum --------------------------------------------------------

# # DEM PREP (run once to generate cudem_3087.tif, then upload to S3)
# # Requires access to T:/05_GIS/TBEP/TBCMP/TBCMP_TBEP_DATASETS.gdb
#
# library(terra)
# gdb <- "T:/05_GIS/TBEP/TBCMP/TBCMP_TBEP_DATASETS.gdb"
# cudem_raw <- rast(paste0('OpenFileGDB:"', gdb, '":NOAA_CUDEM_tbcmp_extent_m'))
#
# # Crop to study area and aggregate to 10m resolution
# load(file = here('data', 'tbcmp_cnt.RData'))
# cudem_raw <- crop(cudem_raw, vect(tbcmp_cnt)) |>
#   mask(vect(tbcmp_cnt)) |>
#   aggregate(fact = 5, fun = "mean")  # 2m -> 10m
#
# # Save locally then upload to S3 for future use
# writeRaster(cudem_raw, here("data/cudem_3087.tif"), filetype = "COG", overwrite = TRUE)
# # aws.s3::put_object(file = here("data/cudem_3087.tif"),
# #                    object = "cudem_3087.tif", bucket = "tbcmp")

# Load county boundaries
load(file = here('data', 'tbcmp_cnt.RData'))

# Load CUDEM topobathy (10m, EPSG:3087, NAVD88 meters) from S3
cudem <- rast("/vsicurl/https://tbcmp.s3.amazonaws.com/cudem_3087.tif")

# Prepare county boundaries
tbcmp_cnt_4326 <- st_transform(tbcmp_cnt, 4326)
bbox_4326 <- st_bbox(tbcmp_cnt_4326)

# Query VDatum API for MLLW -> NAVD88 offsets on a ~5 km sample grid
grid_pts <- build_vdatum_grid(bbox_4326, tbcmp_cnt_4326)
grid_pts <- batch_query_vdatum(grid_pts)

# Interpolate spatially varying MLLW surface aligned to DEM
mllw_surface <- build_mllw_surface(grid_pts, bbox_4326, tbcmp_cnt, cudem)

# Delineate coastal stratum (MLLW to 5 ft NAVD88)
coastal_stratum <- make_coastal_stratum(cudem, mllw_surface)

save(
  coastal_stratum,
  file = here('data', 'coastal_stratum.RData'),
  compress = 'xz'
)

# load files for comparison
load(file = 'T:/04_STAFF/MARCUS/03_GIT/hmpu-workflow/data/coastal.RData')
coastal_4326 <- st_transform(coastal, 4326) |>
  st_zm()
load(file = here('data', 'coastal_stratum.RData'))
coastal_stratum_4326 <- st_transform(coastal_stratum, 4326)

leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    data = coastal_4326,
    color = "blue",
    fillOpacity = 0.5,
    group = "Original TB Stratum"
  ) |>
  addPolygons(
    data = coastal_stratum_4326,
    color = "red",
    fillOpacity = 0.5,
    group = "Interpolated Stratum"
  ) |>
  addLayersControl(
    overlayGroups = c("Original TB Stratum", "Interpolated Stratum"),
    options = layersControlOptions(collapsed = FALSE)
  )

# salinity ---------------------------------------------------------------

# Load county boundaries
load(file = here('data', 'tbcmp_cnt.RData'))

wqp_salinity <- fetch_wqp_salinity(
  counties = tbcmp_cnt,
  by_county = TRUE
)

salinity_layer <- build_salinity_layer(
  sal_pts = wqp_salinity,
  counties = tbcmp_cnt
)

save(
  salinity_layer,
  file = here('data', 'salinity_layer.RData'),
  compress = 'xz'
)

tomap <- wqp_salinity |>
  st_transform(4326)

leaflet(tomap) |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addCircleMarkers(
    color = ~ ifelse(salinity_mean < 5, 'blue', 'red'),
    radius = 4,
    fillOpacity = 0.7,
    label = ~ paste0('Salinity: ', round(salinity_mean, 2))
  )

# map salinity_layer
load(file = here('data', 'salinity_layer.RData'))

salinity_layer_4326 <- st_transform(salinity_layer, 4326)

leaflet(salinity_layer_4326) |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    fillOpacity = 0.7,
    weight = 0.5
  )

# compare to original TB salinity layer
load(file = 'T:/04_STAFF/MARCUS/03_GIT/hmpu-workflow/data/salin.RData')

salin_4326 <- st_transform(salin, 4326)

leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    data = salin_4326,
    fillOpacity = 0.7,
    weight = 0.5,
    color = 'blue',
    group = 'Original TB Salinity'
  ) |>
  addPolygons(
    data = salinity_layer_4326,
    fillOpacity = 0.7,
    weight = 0.5,
    color = 'red',
    group = 'Interpolated Salinity'
  ) |>
  addLayersControl(
    overlayGroups = c('Original TB Salinity', 'Interpolated Salinity'),
    options = layersControlOptions(collapsed = FALSE)
  )

# soils ------------------------------------------------------------------

# Retrieve SSURGO soils data for the 7-county TBCMP area via NRCS Soil Data
# Access (SDA) and classify map units as hydric, mesic, or xeric, matching
# the methods in Ries and Scheda (2014).
#
# Classification:
#   gridcode 100 = Xeric  (well-drained upland soils)
#   gridcode 200 = Mesic  (moderately drained transitional soils)
#   gridcode 300 = Hydric (hydric-rated or poorly drained soils)
#
# Queries are issued county by county to stay within SDA size limits.
# Component attributes are chunked in groups of 500 mukeys.

load(file = here('data', 'tbcmp_cnt.RData'))

soils <- build_soils_layer(tbcmp_cnt, crs = prj)

save(soils, file = here('data', 'soils.RData'), compress = 'xz')

# load original TB-only soils for comparison

load(file = here('data', 'soils.RData'))
soils_4326 <- st_transform(soils, 4326)

soils_tb <- get(load(
  'T:/04_STAFF/MARCUS/03_GIT/hmpu-workflow/data/soils.RData'
)) |>
  st_transform(4326)

tbcmp_cnt_4326 <- st_transform(tbcmp_cnt, 4326)

soil_pal <- colorFactor(
  palette = c('#d4a84b', '#6baed6', '#74c476'),
  levels = c(100, 200, 300)
)

leaflet() |>
  # addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    data = soils_tb,
    color = ~ soil_pal(gridcode),
    fillColor = ~ soil_pal(gridcode),
    fillOpacity = 0.5,
    weight = 0.5,
    label = ~ paste0('TB: ', Descrip),
    group = 'TB Soils (original)'
  ) |>
  addPolygons(
    data = soils_4326,
    color = ~ soil_pal(gridcode),
    fillColor = ~ soil_pal(gridcode),
    fillOpacity = 0.5,
    weight = 0.5,
    label = ~ paste0('TBCMP: ', Descrip),
    group = 'TBCMP Soils (new)'
  ) |>
  addPolygons(
    data = tbcmp_cnt_4326,
    color = 'black',
    fillOpacity = 0,
    weight = 1.5,
    label = ~county,
    group = 'Counties'
  ) |>
  addLegend(
    position = 'bottomright',
    pal = soil_pal,
    values = c(100, 200, 300),
    labFormat = labelFormat(
      transform = function(x) {
        c('Xeric', 'Mesic', 'Hydric')[match(x, c(100, 200, 300))]
      }
    ),
    title = 'Soil Type'
  ) |>
  addLayersControl(
    overlayGroups = c('TB Soils (original)', 'TBCMP Soils (new)', 'Counties'),
    options = layersControlOptions(collapsed = FALSE)
  )

# FNAI -------------------------------------------------------------------

# LULC -------------------------------------------------------------------
