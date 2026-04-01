# setup ------------------------------------------------------------------

library(tidyverse)
library(sf)
library(here)
library(units)
library(flextable)

source(here('R', 'funcs.R'))

# Load shared inputs
load(here('data', '01_inputs', 'fluccs.RData'))
load(here('data', '01_inputs', 'strata.RData'))
load(here('data', '01_inputs', 'tbcmp_cnt.RData'))
load(here('data', '01_inputs', 'coastal_stratum.RData'))

out_dir <- here('data', '03_current_tables')

# current table by county ------------------------------------------------

for (county in tbcmp_cnt$county) {
  county_lower <- tolower(county)
  message('\n--- ', county, ' ---')

  # Load county LULC and seagrass (subtidal)
  load(here('data', '01_inputs', paste0('lulc_', county_lower, '.RData')))
  load(here('data', '01_inputs', paste0('seagrass_', county_lower, '.RData')))
  lulc <- get(paste0('lulc_', county_lower))
  subt <- get(paste0('seagrass_', county_lower))

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
  nativelyr <- get(paste0('nativelyr_', county_lower))
  restorelyr <- get(paste0('restorelyr_', county_lower))

  # Build table
  tab <- curex_fun(
    lulc = lulc,
    subt = subt,
    coastal_stratum = coastal_stratum,
    fluccs = fluccs,
    strata = strata,
    nativelyr = nativelyr,
    restorelyr = restorelyr,
    tbcmp_cnt = tbcmp_cnt,
    county = county,
    cap = paste(
      'Current habitat extent and conservation opportunity -',
      county,
      'County'
    )
  )

  obj_name <- paste0('current_table_', county_lower)
  assign(obj_name, tab)
  save(
    list = obj_name,
    file = file.path(out_dir, paste0(obj_name, '.RData')),
    compress = 'xz'
  )
  message('  Saved as ', obj_name, '.RData')
}
