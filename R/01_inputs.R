# setup ------------------------------------------------------------------

library(tidyverse)
library(sf)
library(curl)
library(here)
library(FedData)

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

# requires DEM
# county boundaries
# VDatum

# salinity ---------------------------------------------------------------

# soils ------------------------------------------------------------------

# FNAI -------------------------------------------------------------------

# LULC -------------------------------------------------------------------
