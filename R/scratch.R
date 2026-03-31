library(tidyverse)
library(sf)
library(here)
library(soilDB)
library(leaflet)

source(here('R/funcs.R'))

prj <- 3087

load(file = '~/Desktop/sg1.RData')
load(file = '~/Desktop/sg2.RData')

leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(
    data = sg1,
    fillOpacity = 0.6,
    color = '#555555',
    weight = 0.5,
    label = ~FLUCCSCODE,
    group = 'sg1'
  ) |>
  addPolygons(
    data = sg2,
    fillOpacity = 0.6,
    color = '#555555',
    weight = 0.5,
    label = ~FLUCCSCODE,
    group = 'sg2'
  ) |>
  addLayersControl(
    overlayGroups = c('sg1', 'sg2'),
    options = layersControlOptions(collapsed = FALSE)
  )
