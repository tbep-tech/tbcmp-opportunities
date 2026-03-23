library(tidyverse)
library(sf)
library(here)
library(soilDB)
library(leaflet)

source(here('R/funcs.R'))

prj <- 3087

# county only -----------------------------------------------------------

load(file = here('data', 'tbcmp_cnt.RData'))

county <- tbcmp_cnt |> filter(county == 'Hillsborough')

# --- 1-4. Fetch, classify, clip, dissolve --------------------------------

soils_county <- build_soils_layer(county, crs = prj)

# --- Compare with original TB soils (clipped to county) ------------------

soils_tb <- get(load(
  'T:/04_STAFF/MARCUS/03_GIT/hmpu-workflow/data/soils.RData'
)) |>
  st_transform(prj) |>
  st_make_valid() |>
  st_intersection(sf::st_union(county)) |>
  st_transform(4326)

soils_county_4326 <- st_transform(soils_county, 4326)

soil_pal <- colorFactor(
  palette = c('#d4a84b', '#6baed6', '#74c476'),
  levels = c(100, 200, 300)
)

leaflet() |>
  # addTiles() |>
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
    data = soils_county_4326,
    color = ~ soil_pal(gridcode),
    fillColor = ~ soil_pal(gridcode),
    fillOpacity = 0.5,
    weight = 0.5,
    label = ~ paste0('TBCMP: ', Descrip),
    group = 'TBCMP Soils (new)'
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
    overlayGroups = c('TB Soils (original)', 'TBCMP Soils (new)'),
    options = layersControlOptions(collapsed = FALSE)
  )
