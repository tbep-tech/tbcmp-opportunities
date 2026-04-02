# setup ------------------------------------------------------------------

library(tidyverse)
library(sf)
library(here)

source(here('R', 'funcs.R'))

# Load shared inputs
load(here('data', '01_inputs', 'tbcmp_cnt.RData'))

out_dir <- here('data', '05_restoration_maps')

# restoration layers by county -------------------------------------------

for (county in tbcmp_cnt$county) {
  county_lower <- tolower(county)
  message('\n--- ', county, ' ---')

  # Load county restorable layer (from 02_current_layers.R)
  load(here(
    'data',
    '02_current_layers',
    paste0('restorelyr_', county_lower, '.RData')
  ))
  restorelyr <- get(paste0('restorelyr_', county_lower))

  obj_name <- paste0('restmap_', county_lower)
  assign(obj_name, restdat_fun(restorelyr = restorelyr))

  # Save as RData
  save(
    list = obj_name,
    file = file.path(out_dir, paste0(obj_name, '.RData')),
    compress = 'xz'
  )

  # Save as shapefile
  sf::st_write(
    get(obj_name),
    file.path(out_dir, paste0(obj_name, '.shp')),
    delete_layer = TRUE
  )

  message('  Saved ', obj_name)
}

# view map
load(here('data', '05_restoration_maps', 'restmap_pinellas.RData'))

restmap_leaflet(restmap_pinellas)
