# Environment preparation ####
rm(list = ls()) # Reset R`s brain

# Load libraries
library(tidyr)
library(dplyr)
library(sf)

# Make a custom operator "not in" for `dplyr` filtering
`%notin%` <- Negate(`%in%`)

path_gbif_sf_dataset <- "./outputs/"

# Vector of datasetKeys for the dataset occurrences from which we 
# deliberately drop from the data
dropped_datasets <- c("c779b049-28f3-4daf-bbf4-0a40830819b6") # EBCC Atlas of European Breeding Birds

# Limit for coordinate uncertainty in meters (occurrences with uncertainty above 
# the threshold will be dropped)
coordUncert.threshold <- 500

# Minimum coordinate precision (occurrences with precision above the threshold will be dropped)
# https://dwc.tdwg.org/terms/#dwc:coordinatePrecision
coordinatePrec.threshold <- 0.001

# load country polygon
country_polygon <- st_read("./shp/country.shp")

# Load data saved at step #1
load(file = "./temp/matches.Rdata")
list2env(matches, .GlobalEnv)
rm(matches)

# Restore original input data, but with internal IDs
attributes <- goodmatch %>% bind_rows(badmatch) %>% 
  arrange(ID)

# Drop unused objects
rm(goodmatch, badmatch)

# Load GBIF data prepared at the second step
load(file = "./temp/gbif_data.Rdata")

# Transform data
gbif_sf_dataset <- gbif.dump %>% 
  # rename columns
  rename(Latitude = decimalLatitude, Longitude = decimalLongitude) %>%
  # make new columns for occurrence and dataset Keys
  mutate(URL_record = paste0("https://www.gbif.org/occurrence/", gbifID)) %>%
  mutate(URL_dataset = paste0("https://www.gbif.org/dataset/", datasetKey)) %>% 
  # drop records from grid datasets
  filter(datasetKey %notin% dropped_datasets) %>% 
  # drop occs w/ high georeference uncertainty
  filter(coordinateUncertaintyInMeters <= coordUncert.threshold | is.na(coordinateUncertaintyInMeters)) %>%
  filter(coordinatePrecision < coordinatePrec.threshold | is.na(coordinatePrecision)) %>% 
  # merge occurrences with attributes (local names, conservation lists)
  left_join(attributes, by = "ID") %>% 
  # drop redundant columns
  select(-c(ID,
            occurrenceID,
            gbifID,
            datasetKey,
            taxonKey.x,
            scientificName.x,
            verbatimScientificName.x,
            taxonKey.y,
            scientificName.y,
            status,
            matchType)) %>% 
  # leave only one scientific name column, with the name as in the input data
  rename(scientificName = verbatimScientificName.y) %>% 
  # coalesce IUCN Red List categories across all sources category) and manually filled (IUCN))
  unite("iucnRListCat", c(iucnRedListCategory, IUCN), sep = "", na.rm = TRUE, remove = FALSE) %>% 
  select(-c(iucnRedListCategory, IUCN)) %>% 
  rename(iucnRedListCategory = iucnRListCat) %>%
  # convert to `simple features` spatial object
  st_as_sf(dim = "XY", remove = FALSE, na.fail = F, 
                      coords = c("Longitude", "Latitude"), crs = "+proj=longlat +datum=WGS84 +no_defs") %>% 
  # clip by polygon and drop unused column
  filter(st_intersects(geometry, country_polygon, sparse = FALSE))

# Preview result
library(ggplot2)
ggplot() +
  geom_sf(data = country_polygon, fill = "transparent", colour = "gray") +
  geom_sf(data = gbif_sf_dataset, colour = "red") +
  theme_bw()

# Save GBIF points to local drive as Robject
save(gbif_sf_dataset, file = "./outputs/gbif_sf_dataset.Rdata")


# Delete temporary files if it exist ####
filestodelete <- list.files(path = "./temp")

for (i in 1:length(filestodelete)) {
  if (file.exists(filestodelete[i])) {
    file.remove(filestodelete[i])
    cat("File deleted")
  } else {
    cat("No file found")
  }
}

# Clean the session
rm(list = ls())
gc()