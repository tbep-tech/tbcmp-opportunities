# setup ------------------------------------------------------------------

library(tidyverse)
library(sf)
library(here)
library(leaflet)

source(here('R', 'funcs.R'))

# Load shared inputs (county-wide layers built in 01_inputs.R)
load(here('data', '01_inputs', 'fluccs.RData'))
load(here('data', '01_inputs', 'tbcmp_cnt.RData'))
load(here('data', '01_inputs', 'coastal_stratum.RData'))
load(here('data', '01_inputs', 'soils.RData'))
load(here('data', '01_inputs', 'salinity_layer.RData'))
load(here('data', '01_inputs', 'prop.RData'))
load(here('data', '01_inputs', 'exst.RData'))

out_dir <- here('data', '02_current_layers')

# current layers by county -----------------------------------------------

# For each county:
#   1. Load county LULC
#   2. Clip shared layers to county boundary
#   3. Build nativelyr, restorelyr, nativersrv, restorersrv
#   4. Save each as <layer>_<county>.RData

for (county in tbcmp_cnt$county) {
  county_lower <- tolower(county)
  message('\n--- ', county, ' ---')

  # Load county LULC
  load(here('data', '01_inputs', paste0('lulc_', county_lower, '.RData')))
  lulc <- get(paste0('lulc_', county_lower))

  # Clip shared layers to county boundary
  cnt_geom <- sf::st_union(tbcmp_cnt[tbcmp_cnt$county == county, ])
  coastal_cnt <- sf::st_intersection(coastal_stratum, cnt_geom)
  soils_cnt <- sf::st_intersection(soils, cnt_geom)
  salin_cnt <- sf::st_intersection(salinity_layer, cnt_geom)
  prop_cnt <- sf::st_intersection(prop, cnt_geom)
  exst_cnt <- sf::st_intersection(exst, cnt_geom)

  # Build the four current-condition layers
  out <- build_current_lyrs(
    lulc = lulc,
    coastal = coastal_cnt,
    soils = soils_cnt,
    salin = salin_cnt,
    prop = prop_cnt,
    exst = exst_cnt,
    fluccs = fluccs
  )

  # Save each layer with county suffix
  for (lyr in names(out)) {
    obj_name <- paste0(lyr, '_', county_lower)
    assign(obj_name, out[[lyr]])
    save(
      list = obj_name,
      file = file.path(out_dir, paste0(obj_name, '.RData')),
      compress = 'xz'
    )
    message('  Saved as ', obj_name, '.RData')
  }
}

# compare Hillsborough new with HMPU
load(file = 'T:/04_STAFF/MARCUS/03_GIT/hmpu-workflow/data/restorelyr.RData')
load(file = here('data', '02_current_layers', 'restorelyr_hillsborough.RData'))
load(file = here('data', '01_inputs', 'tbcmp_cnt.RData'))

cmp1 <- st_transform(restorelyr, 4326)
cmp2 <- st_transform(restorelyr_hillsborough, 4326)

leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    data = cmp1,
    fillOpacity = 0.6,
    color = '#555555',
    weight = 0.5,
    group = 'old'
  ) |>
  addPolygons(
    data = cmp2,
    fillOpacity = 0.6,
    color = '#555555',
    weight = 0.5,
    group = 'new'
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
    overlayGroups = c('old', 'new', 'Counties'),
    options = layersControlOptions(collapsed = FALSE)
  )
