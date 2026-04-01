# setup ------------------------------------------------------------------

library(tidyverse)
library(sf)
library(here)

source(here('R', 'funcs.R'))

# Load shared inputs
load(here('data', '01_inputs', 'tbcmp_cnt.RData'))
load(here('data', '01_inputs', 'coastal_stratum.RData'))

out_dir <- here('data', '04_opportunities_maps')

# opportunity layers by county -------------------------------------------

for (county in tbcmp_cnt$county) {
  county_lower <- tolower(county)
  message('\n--- ', county, ' ---')

  # Load county current layers (from 02_current_layers.R)
  load(here(
    'data',
    '02_current_layers',
    paste0('nativelyr_', county_lower, '.RData')
  ))
  load(here(
    'data',
    '02_current_layers',
    paste0('restorelyr_', county_lower, '.RData')
  ))
  load(here(
    'data',
    '02_current_layers',
    paste0('nativersrv_', county_lower, '.RData')
  ))
  load(here(
    'data',
    '02_current_layers',
    paste0('restorersrv_', county_lower, '.RData')
  ))

  nativelyr <- get(paste0('nativelyr_', county_lower))
  restorelyr <- get(paste0('restorelyr_', county_lower))
  nativersrv <- get(paste0('nativersrv_', county_lower))
  restorersrv <- get(paste0('restorersrv_', county_lower))

  obj_name <- paste0('oppmap_', county_lower)
  assign(
    obj_name,
    oppdat_fun(
      nativersrv = nativersrv,
      restorersrv = restorersrv,
      nativelyr = nativelyr,
      restorelyr = restorelyr,
      coastal_stratum = coastal_stratum,
      tbcmp_cnt = tbcmp_cnt,
      county = county
    )
  )

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
load(here('data', '04_opportunities_maps', 'oppmap_pinellas.RData'))

oppmap_leaflet(oppmap_pinellas)
