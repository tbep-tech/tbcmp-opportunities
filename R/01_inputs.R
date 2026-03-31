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
library(smoothr)
library(rmapshaper)

source(here('R/funcs.R'))

prj <- 3087

# FLUCCS lookup ----------------------------------------------------------

fluccs <- read.csv(
  here('data-raw', '01_inputs', 'FLUCCShabsclass.csv'),
  stringsAsFactors = F
)

save(fluccs, file = here('data', '01_inputs', 'fluccs.RData'), compress = 'xz')

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

save(strata, file = here('data', '01_inputs', 'strata.RData'), compress = 'xz')

# tbcmp counties ---------------------------------------------------------

# # TBCMP project area bounded by 7 county census areas including associated coastal areas tied to local/state jurisdictions
# curl_download(url = "https://www2.census.gov/geo/tiger/TIGER2025/COUSUB/tl_2025_12_cousub.zip",
#              destfile = "./data-raw/01_inputs/tl_2025_12_cousub.zip")
# unzip("./data-raw/01_inputs/tl_2025_12_cousub.zip", exdir = "./data-raw/01_inputs")

fl_counties <- st_read("./data-raw/01_inputs/tl_2025_12_cousub.shp")

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

save(
  tbcmp_cnt,
  file = here('data', '01_inputs', 'tbcmp_cnt.RData'),
  compress = 'xz'
)

# coastal stratum --------------------------------------------------------

# # DEM PREP (run once to generate cudem_3087.tif, then upload to S3)
# # Requires access to T:/05_GIS/TBEP/TBCMP/TBCMP_TBEP_DATASETS.gdb
#
# library(terra)
# gdb <- "T:/05_GIS/TBEP/TBCMP/TBCMP_TBEP_DATASETS.gdb"
# cudem_raw <- rast(paste0('OpenFileGDB:"', gdb, '":NOAA_CUDEM_tbcmp_extent_m'))
#
# # Crop to study area and aggregate to 10m resolution
# load(file = here('data', '01_inputs','tbcmp_cnt.RData'))
# cudem_raw <- crop(cudem_raw, vect(tbcmp_cnt)) |>
#   mask(vect(tbcmp_cnt)) |>
#   aggregate(fact = 5, fun = "mean")  # 2m -> 10m
#
# # Save locally then upload to S3 for future use
# writeRaster(cudem_raw, here("data/01_inputs/cudem_3087.tif"), filetype = "COG", overwrite = TRUE)
# # aws.s3::put_object(file = here("data/01_inputs/cudem_3087.tif"),
# #                    object = "cudem_3087.tif", bucket = "tbcmp")

# Load county boundaries
load(file = here('data', '01_inputs', 'tbcmp_cnt.RData'))

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
  file = here('data', '01_inputs', 'coastal_stratum.RData'),
  compress = 'xz'
)

# load files for comparison
load(file = 'T:/04_STAFF/MARCUS/03_GIT/hmpu-workflow/data/coastal.RData')
coastal_4326 <- st_transform(coastal, 4326) |>
  st_zm()
load(file = here('data', '01_inputs', 'coastal_stratum.RData'))
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
load(file = here('data', '01_inputs', 'tbcmp_cnt.RData'))

salinity_layer <- build_salinity_layer(
  counties = tbcmp_cnt,
  by_county = TRUE
)

save(
  salinity_layer,
  file = here('data', '01_inputs', 'salinity_layer.RData'),
  compress = 'xz'
)

# compare to original TB salinity layer
load(file = 'T:/04_STAFF/MARCUS/03_GIT/hmpu-workflow/data/salin.RData')

salinity_layer_4326 <- st_transform(salinity_layer, 4326)
salin_4326 <- st_transform(salin, 4326)

sal_pal <- colorFactor(
  palette = c('#2166ac', '#74c476', '#d73027'), # blue / green / red
  domain = c('Fresh (<0.5)', '0.5-18', '>18')
)

leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    data = salinity_layer_4326,
    fillColor = ~ sal_pal(Descrip),
    fillOpacity = 0.6,
    color = '#555555',
    weight = 0.5,
    label = ~Descrip,
    group = 'Interpolated Salinity'
  ) |>
  addPolygons(
    data = salin_4326,
    fillOpacity = 0,
    color = '#000000',
    weight = 2.5,
    dashArray = '6 4',
    group = 'Original TB Salinity'
  ) |>
  addLegend(
    pal = sal_pal,
    values = c('Fresh (<0.5)', '0.5-18', '>18'),
    title = 'Salinity zone<br>(Interpolated)',
    position = 'bottomright'
  ) |>
  addLayersControl(
    overlayGroups = c('Interpolated Salinity', 'Original TB Salinity'),
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

load(file = here('data', '01_inputs', 'tbcmp_cnt.RData'))

soils <- build_soils_layer(tbcmp_cnt, crs = prj)

save(soils, file = here('data', '01_inputs', 'soils.RData'), compress = 'xz')

# load original TB-only soils for comparison

load(file = here('data', '01_inputs', 'soils.RData'))
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

# proposed and conservation lands ----------------------------------------

# All FNAI zips are cached to data-raw/01_inputs/fnai/ so re-runs skip the download.
# Update the filenames in the URLs when FNAI publishes newer versions.

load(file = here('data', '01_inputs', 'tbcmp_cnt.RData'))

# FLMA: Florida Conservation Lands (excludes MacDill AFB)
flma <- fetch_fnai(
  url = 'https://www.fnai.org/shapefiles/flma_202503.zip',
  cnt = tbcmp_cnt
) |>
  dplyr::filter(!MANAME %in% 'MacDill Air Force Base')

save(flma, file = here('data', '01_inputs', 'flma.RData'), compress = 'xz')

# FFBOT: Future Forever Board of Trustees projects
ffbot <- fetch_fnai(
  url = 'https://www.fnai.org/shapefiles/ffbot_202503.zip',
  cnt = tbcmp_cnt
)

save(ffbot, file = here('data', '01_inputs', 'ffbot.RData'), compress = 'xz')

# FFA: Future Forever acquisitions
ffa <- fetch_fnai(
  url = 'https://www.fnai.org/shapefiles/ff_acquired_202502.zip',
  cnt = tbcmp_cnt
)

save(ffa, file = here('data', '01_inputs', 'ffa.RData'), compress = 'xz')

# Aquatic Preserves: Florida DEP
# Source: https://geodata.dep.state.fl.us/datasets/81841412d3984e9aac2c00c21e41d32e_0

aqprs <- fetch_aqprs(cnt = tbcmp_cnt)

save(aqprs, file = here('data', '01_inputs', 'aqprs.RData'), compress = 'xz')

# CLIP: Critical Lands and Waters Identification Project (FNAI)
# Source: https://www.fnai.org/services/clip
# Raster inside File GDB read via GDAL OpenFileGDB driver (requires GDAL >= 3.7).
# Values 1-5; 1 = highest priority. Zip removed after processing.

load(file = here('data', '01_inputs', 'tbcmp_cnt.RData'))

clip <- fetch_clip(cnt = tbcmp_cnt)

save(clip, file = here('data', '01_inputs', 'clip.RData'), compress = 'xz')

# compare to existing TB CLIP layer
load(file = here('data', '01_inputs', 'clip.RData'))
clip_old <- sf::st_read(
  'https://opendata.arcgis.com/datasets/ba464bae7a1144f09522b459d297d1ef_10.geojson'
) |>
  sf::st_transform(4326) |>
  st_make_valid() |>
  dplyr::group_by(gridcode) |>
  dplyr::summarise(geometry = sf::st_union(geometry), .groups = 'drop') |>
  dplyr::rename(priority = gridcode)

clip_4326 <- clip |>
  dplyr::filter(priority >= 1) |>
  sf::st_transform(4326)

clip_pal <- colorFactor(
  palette = c('#006837', '#78c679', '#d9f0a3'),
  domain = 1:3
)

leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    data = clip_4326,
    fillColor = ~ clip_pal(priority),
    fillOpacity = 0.6,
    color = '#555555',
    weight = 0.5,
    label = ~ paste('Priority', priority),
    group = 'CLIP v4 (TBCMP)'
  ) |>
  addPolygons(
    data = clip_old,
    fillOpacity = 0,
    color = '#000000',
    weight = 2,
    dashArray = '6 4',
    group = 'Existing TB CLIP'
  ) |>
  addLegend(
    pal = clip_pal,
    values = 1:3,
    title = 'CLIP Priority',
    position = 'bottomright',
    labFormat = labelFormat(prefix = 'Priority ')
  ) |>
  addLayersControl(
    overlayGroups = c('CLIP v4 (TBCMP)', 'Existing TB CLIP'),
    options = layersControlOptions(collapsed = FALSE)
  )

# proposed and existing conservation lands --------------------------------

load(file = here('data', '01_inputs', 'flma.RData'))
load(file = here('data', '01_inputs', 'ffbot.RData'))
load(file = here('data', '01_inputs', 'ffa.RData'))
load(file = here('data', '01_inputs', 'aqprs.RData'))
load(file = here('data', '01_inputs', 'clip.RData'))

out <- build_prop_exst(flma, ffbot, ffa, aqprs, clip)
prop <- out$prop
exst <- out$exst

save(prop, file = here('data', '01_inputs', 'prop.RData'), compress = 'xz')
save(exst, file = here('data', '01_inputs', 'exst.RData'), compress = 'xz')

# LULC -------------------------------------------------------------------

# save 2023 lulc files to data folder by county
fetch_lulc()

# check county lulc
load(file = here('data', '01_inputs', 'lulc_pinellas.RData'))
load(file = here('data', '01_inputs', 'lulc_hillsborough.RData'))

lulc_pinellas_4326 <- st_transform(lulc_pinellas, 4326)
lulc_hillsborough_4326 <- st_transform(lulc_hillsborough, 4326)

leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    data = lulc_pinellas_4326,
    fillOpacity = 0.6,
    color = '#555555',
    weight = 0.5,
    label = ~FLUCCSCODE,
    group = 'pinellas'
  ) |>
  addPolygons(
    data = lulc_hillsborough_4326,
    fillOpacity = 0.6,
    color = '#555555',
    weight = 0.5,
    label = ~FLUCCSCODE,
    group = 'hillsborough'
  ) |>
  addLayersControl(
    overlayGroups = c('pinellas', 'hillsborough'),
    options = layersControlOptions(collapsed = FALSE)
  )

# seagrass ---------------------------------------------------------------

# Save combined seagrass clipped by county
fetch_seagrass()

# check county seagrass
load(file = here('data', '01_inputs', 'seagrass_manatee.RData'))
load(file = here('data', '01_inputs', 'seagrass_sarasota.RData'))
load(file = here('data', '01_inputs', 'tbcmp_cnt.RData'))

co1_4326 <- st_transform(seagrass_manatee, 4326)
co2_4326 <- st_transform(seagrass_sarasota, 4326)

leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    data = co1_4326,
    fillOpacity = 0.6,
    color = '#555555',
    weight = 0.5,
    label = ~FLUCCSCODE,
    group = 'co1'
  ) |>
  addPolygons(
    data = co2_4326,
    fillOpacity = 0.6,
    color = '#555555',
    weight = 0.5,
    label = ~FLUCCSCODE,
    group = 'co2'
  ) |>
  addPolygons(
    data = st_transform(tbcmp_cnt, 4326),
    fillOpacity = 0,
    color = 'black',
    weight = 1.5,
    label = ~county,
    group = 'Counties'
  ) |>
  addLayersControl(
    overlayGroups = c('co1', 'co2', 'Counties'),
    options = layersControlOptions(collapsed = FALSE)
  )
