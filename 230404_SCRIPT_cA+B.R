# Chapter I: Conservation proposal for the Espeletiinae subtribe (Asteraceae) based on the Important Plant Areas (IPAs) methodology
# IPA methodology (Diazgranados & Castellanos-Castro 2018)

rm(list = ls())

## Installation Packages
packages_list<-list("data.table", "sp", "sf", "fasterize","raster", "rgdal", "stars", "foreign", "terra", "cli", "vegan", "tidyverse", "scales", "plyr", "fdth", "reshape", "reshape2","ggnewscale","raster", "rgeos", "ggplot2", "tmap", "pryr","base", "maptools", "ggmap", "rnaturalearth", "rnaturalearthdata", "rstudioapi", "magrittr", "pbapply", "rredlist", "grid", "gdalUtilities", "rCAT")
lapply(packages_list, library, character.only = TRUE)

# With the coordinates of the records, the study area is delimited, cutting with a DEM layer of 90m from 1000m
# With the section of the study area, a grid of 10 km is established
setwd("C:/Users/Clari Mora/Dropbox/PC/Documents")
dirfile<- "C:/RESULTADOS IPAs_2022/Poligono_Area Estudio_1000 m_UA 10x10km/UA_1000_10x10.shp"

# Call the study area and transform it into a raster with flat coordinates
spatial_file<- st_read(dirfile) %>% st_transform(3395)

# Once loaded, the polygon (working base raster) is rasterized. 
# Creation of identification id for each 10 x 10 km2 pixel
rasterbase<- raster(extent(spatial_file),crs= paste0("+init=epsg:", 3395), res= 10000 )
tname3<- tempfile(fileext = ".tif")
t_file<- writeStart(rasterbase, filename = tname3,  overwrite=T); writeStop(t_file);
gdalUtilities::gdal_rasterize(dirfile, tname3, burn=1, at=T)

# Creating a study area mask
raster_file<- raster(tname3) %>% {terra::mask(setValues(., seq(ncell(.))), .)} %>% rast

# Transformation of the mask into coordinate system 4326
area_4326<- raster(raster_file) %>% projectRaster(crs = st_crs(4326)$proj4string, method = "ngb")
cell_area<- terra::cells(raster_file)

# Generation of study area polygon
poly =rasterToPolygons(raster(raster_file)) %>% st_as_sf() %>% st_transform (4326) 
names(poly) [1]="Id"
st_write(poly, "poligono_pixel_all_v4.shp")

## Loading the records of Espeletiinae
glob <- read.csv("C:/RESULTADOS IPAs_2022/221103_OCC_SRTM_Copernicus_TNC_Useful_Coord.csv", sep=",", header = TRUE, row.names = NULL, na.strings = "")

# Crossing the polygon of the study area with registers. Assignment of the pixel id, eliminating the information that does not have NAs 
glob <- glob %>% 
  mutate(LATITUDE= as.numeric(gsub(",", ".", LATITUDE)), LONGITUD= as.numeric(gsub(",", ".", LONGITUD)), id_coord= group_indices_(., .dots=c("LATITUDE", "LONGITUD")) ) %>%
  dplyr::filter(!is.na(LATITUDE) | !is.na(LONGITUD))

# Relate records and study area polygon, assign it to coordinate system
glob_spatial<- glob %>% 
  st_as_sf(coords = c("LONGITUD","LATITUDE"), crs = 4326)

# Display of records in study area polygon
plot(glob_spatial[, "geometry"])

# Generate a spatial file with the study area coordinate system. "Id" is assigned. It is converted to a data.frame to be used as a table to allow the unions and comparisons of the methodology analysis.
glob_spatial_mask<- as.data.frame(glob_spatial) %>% mutate(Id= raster::extract(area_4326, glob_spatial) ) %>%
  dplyr::filter(!is.na(Id))  %>%  st_as_sf() %>% st_transform(3395)

# Data.frame creation for information management
glob_data<- as.data.frame(glob_spatial_mask)

# Perform record counts by species and by different coordinate
# n_occ_sp: # species records per pixel
# n_occ_sp_pixel: # sp records per pixel. sp count, per pixel, and amount per pixel (sum of records/pixel)
data_id_sp<- glob_data %>% group_by(Id,name_clean, PAIS) %>% dplyr::summarise(n_occ_sp_pixel = n_distinct(id_coord))

# Glob data grouped by sp name, how many times is the sp in total and per pixel joined by name_clean
# sp by Id/pixel, and occ such sp / Id. # occ of the sp: how many sp are there in total
data_sp<- data_id_sp %>% group_by(name_clean) %>% dplyr::summarise(n_occ_sp_total = sum(n_occ_sp_pixel))

# List of total records = records, number of species per Id/pixel [analysis basis].
data_summ<- list(data_id_sp, data_sp) %>% join_all()

# Export result
write.table(data_summ,file="occ_sp_pixel_v2.csv", sep = ",", row.names = TRUE, col.names=TRUE)

##############################################################################
## CRITERIO cA1: sp GLOBALLY THREATENED
# The area contains important populations of globally threatened taxa
# Condition: at least 1% of all recorded specimens per globally threatened species are within a pixel
##############################################################################

# Filter by threatened column, under IUCN columns. Create unique vector with required values (threat sp)
glob_sp_IUCN<- dplyr::filter(glob_data, IUCN_RL_CO %in% c("VU", "EN", "CR"))$name_clean %>% unique()

# Threshold counts per species over total pixels occupied
# How many occurrences there are per sp in the whole mask and then carlcular that % (1%) represents and if classified according to the threshold
glob_data_IUCN_1perc<- dplyr::filter(data_summ, name_clean %in% glob_sp_IUCN) %>% 
  mutate(percent= n_occ_sp_pixel/ n_occ_sp_total) %>% ## add nueva columna con el calculo con el % calculado del umbral[mutate]
  dplyr::filter(percent>=0.01)

# by a gropup_by, per pixel Id, how many sp with at least 1% meet
glob_data_IUCN_1perc_spatial <- glob_data_IUCN_1perc %>% group_by(Id) %>% 
  dplyr::summarise(cA1= n_distinct(name_clean), name_clean=name_clean) %>% 
  distinct() %>% list(dplyr::select(data_summ, -n_occ_sp_total)) %>% join_all

# Creating table with Id pixel, cA1 and occ per pixel columns
glob_data_IUCN_1perc_spatial_uniques <- dplyr::select(glob_data_IUCN_1perc_spatial, -name_clean) %>% distinct()

# Mapping of result and export
rastercA1<- rasterbase
rastercA1[glob_data_IUCN_1perc_spatial_uniques$Id]<- glob_data_IUCN_1perc_spatial_uniques$cA1
plot(raster_file, col= "gray")
plot(rastercA1, add=T)
names(rastercA1) = "cA1"

# Organization of information, by pixel and within this how many sp and how many times are found
# Column cA1_v1: amount of sp that are on each pixel
summ_sp_pixel_cA1 <- split(glob_data_IUCN_1perc_spatial, glob_data_IUCN_1perc_spatial$Id) %>% 
  lapply(function(x) {dplyr::mutate( x[1,], name_clean= paste(x$name_clean, collapse= ";" ) )}) %>%
  rbind.fill()

# Consolidate polygon (shp) and csv base with information - cA1 result
poly_cA1 =rasterToPolygons(rastercA1) %>% st_as_sf() %>% 
  as.data.frame() %>% dplyr::mutate(Id = terra::cells(rast(rastercA1))) %>%
  list(summ_sp_pixel_cA1) %>% join_all() %>% st_as_sf %>%
  st_transform(4326)

st_write(poly_cA1, "cA1_evaluado.shp")
write.table(glob_data_IUCN_1perc_spatial,file="cA1_evaluado_vf.csv", sep = ",", row.names = TRUE, col.names=TRUE)

# data extraction based on counting columns. Display ID pixel column per sp
draft_table_cA1 <- dcast(glob_data_IUCN_1perc_spatial, Id ~ name_clean, drop=TRUE, fill = 0, value.var = "n_occ_sp_pixel")
names(draft_table_cA1)
write.table(draft_table_cA1, file="ccA1_organice.csv", sep = ",", row.names = TRUE, col.names=TRUE)

# Organization and calculation of table percentages by Id and occ by species, for each pixel. 
# Rank: corresponds to the organization of the ranking column, standardized from carried to scale of 0 to 1
percent_IUCN_1perc <- glob_data_IUCN_1perc %>% group_by(Id,name_clean, PAIS) %>%
  dplyr::summarise(perc_occ = n_occ_sp_pixel/ n_occ_sp_total, threshold = n_occ_sp_total * 0.01) %>%
  mutate(cumple = ifelse (perc_occ >= 0.01, "Yes", "No")) %>%
  split(. $name_clean) %>% 
  lapply(function(x) mutate (x, RankcA1 = group_indices_(x, .dots=c("perc_occ"))) %>% mutate (RankcA1 = max(. $RankcA1) -RankcA1+1) %>% 
           mutate(NormcA1 = ((RankcA1 *100) / sum (. $RankcA1)) / 100)
           ) %>% rbind.fill()

write.table(percent_IUCN_1perc, file="cA1_1perc_vf.csv", sep = ",", row.names = TRUE, col.names=TRUE)

##############################################################################
## CRITERION cA3: sp NATIONALLY THREATENED
# The area contains significant populations of nationally threatened taxa
# Condition: at least 10% of all the specimens registered by threatened species in the country, are within a pixel
##############################################################################

# Filter by column of threatened, by country by IUCN value. Create unique vector with required values

country_sp_IUCN <- split(glob_data, glob_data$PAIS) %>% 
  lapply(function(x) dplyr::filter(x, GEPC_Natio %in% c("VU", "EN", "CR"))$name_clean) %>%
  unlist() %>% unique()

# know which are the sp that of the global, are not in the national (it is informative)
sp_global_no_country <- glob_sp_IUCN %>% {.[!. %in% country_sp_IUCN]}
sp_coutry_no_global <- country_sp_IUCN %>% {.[!. %in% glob_sp_IUCN]}

# Counts under Thresholds by species over the total number of occupied pixels
# How many occ are there per sp in the whole mask and then calculate what % (10%) represents and if they classify according to the threshold
country_data_IUCN_10perc <- dplyr::filter(data_summ, name_clean %in% country_sp_IUCN) %>% 
  mutate(percent= n_occ_sp_pixel/ n_occ_sp_total) %>%
  dplyr::filter(percent>=0.1)

# by means of a groupup_by, by the Id columns, how many sp with at least 1% meet to generate a distribution map
country_data_IUCN_10perc_spatial <- country_data_IUCN_10perc %>% group_by(Id) %>% 
  dplyr::summarise(cA3= n_distinct(name_clean), name_clean=name_clean) %>% 
  distinct() %>% list(dplyr::select(data_summ, - n_occ_sp_total)) %>% join_all

# Creation of a table with columns Id pixel, Criteria cA3 and occ per pixel.
country_data_IUCN_10perc_spatial_uniques <- dplyr::select(country_data_IUCN_10perc_spatial, -name_clean) %>% distinct()

# Result mapping and export
rastercA3<- rasterbase
rastercA3[country_data_IUCN_10perc_spatial_uniques$Id]<- country_data_IUCN_10perc_spatial_uniques$cA3
plot (raster_file, col= "gray")
plot(rastercA3, add=T)
names(rastercA3) = "cA3"

# organization information, per pixel and within it how many sp and how many times they are found.
# Column cA3: amount of sp that fall in pixel.
summ_sp_pixel_cA3 <- split(country_data_IUCN_10perc_spatial, country_data_IUCN_10perc_spatial$Id) %>% 
  lapply(function(x) {dplyr::mutate( x[1,], name_clean= paste(x$name_clean, collapse= ";" ) )})  %>%
  rbind.fill()

# Consolidate polygon (shp) and base in csv with information - cA3 result
poly_cA3 =rasterToPolygons(rastercA3) %>% st_as_sf() %>% 
  as.data.frame() %>% dplyr::mutate(Id = terra::cells(rast(rastercA3))) %>%
  list(summ_sp_pixel_cA3) %>% join_all() %>% st_as_sf %>%
  st_transform(4326)

# Export
st_write(poly_cA3, "cA3_evaluado.shp")
write.table(country_data_IUCN_10perc_spatial,file="ccA3_evaluado.csv", sep = ",", row.names = TRUE, col.names=TRUE)

# column-based data extraction for analysis
draft_table_cA3 <- dcast(country_data_IUCN_10perc_spatial, Id ~ name_clean, drop=TRUE, fill = 0, value.var = "n_occ_sp_pixel")
write.table(draft_table_cA3, file="ccA3_organice.csv", sep = ",", row.names = TRUE, col.names=TRUE)

# Organization and calculation of table percentages organized by Id and occ by species, for each sp in each pixel
# Percentage, organization of the % (high and low: Sort by %), data normalization (scale 0 to 1)
percent_country_IUCN_10perc <- country_data_IUCN_10perc %>% group_by(Id,name_clean, PAIS) %>%
  dplyr::summarise(perc_occ_country = n_occ_sp_pixel/ n_occ_sp_total, threshold = n_occ_sp_total * 0.1) %>%
  mutate(cumple = ifelse (perc_occ_country >= 0.1, "Yes", "No")) %>%
  split(. $name_clean) %>% 
  lapply(function(x) mutate (x, RankcA3 = group_indices_(x, .dots=c("perc_occ_country"))) %>% mutate (RankcA3 = max(. $RankcA3) -RankcA3+1) %>% 
           mutate(NormcA3 = ((RankcA3 *100) / sum (. $RankcA3)) / 100)
  ) %>% rbind.fill()

write.table(percent_country_IUCN_10perc, file="cA3_10perc_vf.csv", sep = ",", row.names = TRUE, col.names=TRUE)

##############################################################################
##### IUCN EVALUATION - EOO (extent of presence) and AOO (area of occupancy)
#############################################################################

# DB initial input [glob <- read.csv(system.file("extdata", "221103_OCC_SRTM_Copernicus_TNC_Useful_Coord.csv")] filter by column of sp and coord
occ <-data.frame(taxa=glob$name_clean, lat=glob$LATITUDE,long=glob$LONGITUD)


# Result of EOO and AOO calculation, IUCN category assignment
resultsdfCat_IUCN <- ConBatch(occ,cellsize = 2000,"km", FALSE)

write.table(resultsdfCat_IUCN, file="results_IUCN.csv", sep = ",", row.names = TRUE, col.names=TRUE)

##############################################################################
##### CRITERIO cA4: occ EAR
# The area contains important populations of endemic (or near-endemic) taxa of highly restricted range (EAR)
# Conditions: (i) An EAR species is one whose extent of occurrence is equal to or less than 100 km2; (ii) at least 10% of each species in EAR range is found in a pixel
##############################################################################

EOO_EAR = resultsdfCat_IUCN %>% dplyr::mutate(EOOkm2 = as.numeric(EOOkm2))
sp_EAR <- EOO_EAR %>% dplyr::filter (EOOkm2 <= 100)
View(sp_EAR)

# Calculate how many there are per pixel and filter by column of sp. and with the condition of EAR
sp_EAR_summ <- dplyr::filter(data_summ, name_clean %in% sp_EAR$taxa)
write.table(sp_EAR_summ,file="sp_EAR_cA4.csv", sep = ",", row.names = TRUE, col.names=TRUE)

# Threshold calculation. Extract relative frequencies to the lesser 10% of all the occ of sp
sp_EAR_10perc <- sp_EAR_summ %>% group_by(Id) %>% 
  mutate(percent = n_occ_sp_pixel / n_occ_sp_total) %>% 
  dplyr::filter(percent >= 0.1)
view(sp_EAR_10perc)

# pixel, # occ, % representing. Add an extra filter: cA4
sp_EAR_10perc_pixel <- sp_EAR_10perc %>% list(sp_EAR_summ) %>% join_all %>% group_by(Id) %>% 
  dplyr:: mutate(cA4 = n_distinct(name_clean), name_clean = name_clean) %>% 
  distinct() %>% list (dplyr::select(data_summ, - n_occ_sp_total)) %>% join_all
view(sp_EAR_10perc_pixel)

sp_EAR_10perc_unique <- dplyr:: select(sp_EAR_10perc_pixel, -name_clean) %>%  distinct ()
sp_EAR_10perc_unique

rastercA4 <- rasterbase
rastercA4 [sp_EAR_10perc_unique$Id] <- sp_EAR_10perc_unique$cA4
plot(raster_file, col= "gray")
plot(rastercA4, add=T)
names(rastercA4) = "cA4"

summ_sp_pixel_cA4 <- split(sp_EAR_10perc_pixel, sp_EAR_10perc_pixel$Id) %>% 
  lapply(function(x) {dplyr::mutate( x [1,], name_clean= paste(x$name_clean, collapse = ";") )}) %>% 
  rbind.fill()

poly_cA4 = rasterToPolygons(rastercA4) %>% st_as_sf() %>% 
  as.data.frame() %>% dplyr::mutate(Id = terra::cells(rast(rastercA4))) %>% 
  list(summ_sp_pixel_cA4) %>% join_all() %>%  st_as_sf %>% 
  st_transform(4326)

st_write(poly_cA4, "cA4_evaluado.shp")

write.table(sp_EAR_10perc_pix, file="cA4_evaluado_EAR.csv", sep =",", row.names = TRUE, col.names=TRUE)

draft_table_cA4 <- dcast(sp_EAR_10perc_pixel, Id ~ name_clean, drop=TRUE, fill = 0, value.var = "n_occ_sp_pixel")

write.table(draft_table_cA4, file="cA4_EAR_organice.csv", sep = ",", row.names = TRUE, col.names = TRUE)


# Organization and calculation of table percentages organized by Id and occ by species, for each sp in each pixel
# Percentage, organization of the % (high and low: Sort by %), data normalization (scale 0 to 1)
percent_cA4_EAR <- sp_EAR_10perc_pixel %>% group_by(Id,name_clean, PAIS) %>%
  dplyr::summarise(perc_occ_EAR = n_occ_sp_pixel/ n_occ_sp_total, threshold = n_occ_sp_total * 0.1) %>%
  mutate(cumple = ifelse (perc_occ_EAR >= 0.1, "Yes", "No")) %>%
  split(. $name_clean) %>% 
  lapply(function(x) mutate (x, RankcA4 = group_indices_(x, .dots=c("perc_occ_EAR"))) %>% mutate (RankcA4 = max(. $RankcA4) -RankcA4+1) %>% 
           mutate(NormcA4 = ((RankcA4 *100) / sum (. $RankcA4)) / 100)
  ) %>% rbind.fill()

write.table(percent_cA4_EAR, file="cA4_EAR_result_vf.csv", sep = ",", row.names = TRUE, col.names=TRUE)


##############################################################################
##### CRITERIO cA5: occ ERR
# The area contains important populations of endemic (or near-endemic) taxa of restricted range (ERR)
# Conditions: (i) An ERR species is one whose extent of occurrence is greater than 100 km2 but equal to or less than 5000 km2; (ii) at least 10% of each species in ERR range is found in a pixel
##############################################################################

EOO_ERR = resultsdfCat_IUCN %>% dplyr::mutate ( ( (EOOkm2 > 100) & (EOOkm2 <= 5000) ) )
sp_ERR <- EOO_EAR %>% dplyr::filter (EOOkm2 <= 100)
View(sp_ERR)

# Calculate how many there are per pixel and filter by column of sp. and with the condition of EAR
sp_ERR_summ <- dplyr::filter(data_summ, name_clean %in% sp_ERR$taxa)

# n_ERR = unique(sp_ERR_data_summ$name_clean) %>% length()
sp_ERR_occ = unique(sp_ERR$name_clean) %>% length()

# Threshold calculation. Extract relative frequencies to the lesser 10% of all the occ of sp
sp_ERR_10perc <- sp_ERR_summ %>% group_by(Id) %>% 
  mutate(percent= n_occ_sp_pixel/n_occ_sp_total) %>% 
  dplyr::filter(percent >= 0.1)
view(sp_ERR_10perc)

# pixel, # occ, % representing. Add an extra filter: cA5
sp_ERR_10perc_pixel <- sp_ERR_10perc %>% list(sp_ERR_summ) %>% join_all %>% group_by(Id) %>%
  dplyr::mutate(cA5= n_distinct(name_clean), name_clean= name_clean) %>% 
  distinct() %>% list(dplyr::select(data_summ, - n_occ_sp_total)) %>% join_all
View(sp_ERR_10perc_pixel)

sp_ERR_10perc_unique <- dplyr::select(sp_ERR_10perc_pixel, -name_clean) %>% distinct()
sp_ERR_10perc_unique

rastercA5 <- rasterbase
rastercA5 [sp_ERR_10perc_unique$Id] <- sp_ERR_10perc_unique$cA5
plot(raster_file, col= "gray")
plot(rastercA5, add=T)
names(rastercA5) = "cA5"

summ_sp_pixel_cA5 <- split(sp_ERR_10perc_pixel, sp_ERR_10perc_pixel$Id) %>% 
  lapply(function(x) {dplyr::mutate( x[1,], name_clean= paste(x$name_clean, collapse= ";" ) )})  %>%
  rbind.fill()

poly_cA5 =rasterToPolygons(rastercA5) %>% st_as_sf() %>% 
  as.data.frame() %>% dplyr::mutate(Id = terra::cells(rast(rastercA5))) %>%
  list(summ_sp_pixel_cA5) %>% join_all() %>% st_as_sf %>%
  st_transform(4326)

st_write(poly_cA5, "cA5_evaluado.shp")

write.table(sp_ERR_10perc_pix,file="cA5_evaluado_ERR.csv", sep = ",", row.names = TRUE, col.names=TRUE)

draft_table_cA5 <- dcast(sp_ERR_10perc_pixel, Id ~ name_clean, drop=TRUE, fill = 0, value.var = "n_occ_sp_pixel")
write.table(draft_table_cA5, file="ccA5_ERR_organice.csv", sep = ",", row.names = TRUE, col.names=TRUE)

# Organization and calculation of table percentages organized by Id and occ by species, for each sp in each pixel
# Percentage, organization of the % (high and low: Sort by %), data normalization (scale 0 to 1)
percent_cA5_ERR <- sp_ERR_10perc_pixel %>% group_by(Id,name_clean, PAIS) %>%
  dplyr::summarise(perc_occ_ERR = n_occ_sp_pixel/ n_occ_sp_total, threshold = n_occ_sp_total * 0.1) %>%
  mutate(cumple = ifelse (perc_occ_ERR >= 0.1, "Yes", "No")) %>%
  split(. $name_clean) %>% 
  lapply(function(x) mutate (x, RankcA5 = group_indices_(x, .dots=c("perc_occ_ERR"))) %>% mutate (RankcA5 = max(. $RankcA5) -RankcA5+1) %>% 
           mutate(NormcA5 = ((RankcA5 *100) / sum (. $RankcA5)) / 100)
  ) %>% rbind.fill()

write.table(percent_cA5_ERR, file="cA5_ERR_result.csv", sep = ",", row.names = TRUE, col.names=TRUE)

########################
# CRITERIO B
# B1: The area has the highest estimated richness of native plants and fungi in that type of ecosystem
#######################

# create temporary file
tname_ecosystm <- tempfile(fileext = ".tif")

# previously created rasterbase is used with the pixel ids
t_file<- writeStart(rasterbase, filename = tname_ecosystm,  overwrite=T); writeStop(t_file);

# Import ecosystem layer: Sayre, R., 2022, World Terrestrial Ecosystems (WTE) 2020
# U.S. Geological Survey data release, https://doi.org/10.5066/P9DO61LP.io
dir_ecosystems<- "C:/Users/Clari Mora/Dropbox/PC/Desktop/202210_BackUp_Info Clara/Doctorado/ESPELETIA/ArcGIS Insumos/230201_World Ecosystem SHP/World_Ecosystem.shp"

# rasterize
gdalUtilities::gdal_rasterize(dir_ecosystems, tname_ecosystm,  at=T, a= "gridcode")

# replace shp with dbd
library(foreign)
data_ecosystems<- read.dbf(gsub(".shp", ".dbf", dir_ecosystems))[, c("gridcode", "W_Ecosystm")] %>% distinct()

# layer rasterization
rast_ecosystems<- rast(tname_ecosystm)

# convert layer to table (reduce computational expense)
# rast ecosystem: raster layer
# terra:cell: tell which cells have values - raster base package, terra second version of raster
# you must put gridcode: so that it does the union
# which cell corresponds to each grid code (of my pixel system)
# call the two columns of interest from my grid system and ecosystems. I give pixels the value of gridcode. id never change
table_ecosystems <- as.data.frame(rast_ecosystems) %>% setNames("gridcode") %>% mutate(Id= terra::cells(rast_ecosystems)) %>%
  list(data_ecosystems, data_summ) %>% join_all()

# remove the NAs of those sp that do not fall into any ecosystem
# sp per ecosystem per UA
sp_ecosystems <-table_ecosystems %>% dplyr::filter(!is.na(name_clean)) %>% 
  group_by(W_Ecosystm) %>% dplyr::summarise(nsp_ecosystem= n_distinct(name_clean))

# calculation of how many occ there are by ecosystems: how many different sp there are
sp_pixel <- table_ecosystems %>% group_by(Id) %>% dplyr::summarise(nsp_pixel= n_distinct(name_clean), W_Ecosystm=W_Ecosystm)

# ecosystem listing, # sp per pixel and per ecosystem: this is to join sp and ecosystem information
# we have real and projected richness and all those ecosystems that have more than 50% that are reprojected: they are accepted (removing 3, 2 1, records)
# don't filter by a specific number, because that would be arbitrary
sp_pixel_ecosystem <- list(sp_ecosystems, sp_pixel) %>% join_all() %>% mutate(perc= nsp_pixel/ nsp_ecosystem) %>%
  dplyr::filter()

write.table(sp_pixel_ecosystem,file="sp_pixel_W_Ecosystm_v1.csv", sep = ",", row.names = TRUE, col.names=TRUE)

# partition species listings by ecosystem for those with sp
list_ecosystems <- dplyr::filter(table_ecosystems, W_Ecosystm %in% unique(sp_ecosystems$W_Ecosystm) ) %>% 
  dplyr::mutate(W_Ecosystm= factor(W_Ecosystm)) %>% 
  split(.$W_Ecosystm)

# Calculation of Chao (estimated richness) by ecosystem.
chao_ecosystems<- lapply(list_ecosystems, function(x){ print(unique(x$W_Ecosystm))
  
  tryCatch({
    
    draft_ecosystems <- dcast(x, Id ~ name_clean, drop=TRUE, fill = 0, value.var = "n_occ_sp_pixel")  %>%
      {.[, !names(.) %in% "NA"]} %>% column_to_rownames("Id")
    
    rich_ecosystems_estimated <- specpool(draft_ecosystems) %>% mutate(W_Ecosystm= unique(x$W_Ecosystm))
  }, error= function(e) {data.frame(W_Ecosystm=unique(x$W_Ecosystm))})
  
}) %>% rbind.fill() %>% mutate(repres= Species/chao) %>% dplyr::filter(Species>1)




write.csv(chao_ecosystems, file="sp_chao_ecosystems_v1.csv", sep = ",", row.names = TRUE, col.names=TRUE)   

###################
## CRITERION B2
# B2: That area/s with an estimated richness of plants complementary to B1, so that together they add up to 10% of the total biodiversity of the ecosystem
# Condition cB2: when the AIP selected by criterion B1 contains 10% or more of the total richness estimated for the type of ecosystem, the application of criterion B2 will not be necessary.
##################

# column "repres": how representative is the pixel with respect to the total: # sp / in the entire ecosystem (Chao1)
# ceiling: to round decimals, otherwise give full units
# everything with more than 10% will be cB2
chao_richness<- data.frame(chao_ecosystem= chao_ecosystems$chao, W_Ecosystm= chao_ecosystems$W_Ecosystm) %>%
  list(table_ecosystems) %>% join_all() %>% group_by(Id) %>%
  dplyr::filter(!is.na(name_clean)) %>% 
  dplyr::mutate(sp_pixel= n_distinct(name_clean)) %>% dplyr::rowwise() %>%
  dplyr::mutate(chao_ecosystem= ceiling(chao_ecosystem)) %>% dplyr::mutate(repres=sp_pixel / (chao_ecosystem)) %>% 
  dplyr::mutate(B2= ifelse(repres>=0.1, "SI", NA))


# pixels with representativeness > 10%
# display table with bye, Id, name_clean, occ, sp_pixel, representativeness according to bye, and column if it complies with B2
B2_pixel<- dplyr::filter(chao_richness, !is.na(B2)) %>%  data.frame()
write.csv(B2_pixel, file = "cB2_complementario.csv")


# cB2 map output
# starting from an estimated richness per ecosystem, it is indicated that these are the most important pixels per ecosystem due to the data they have
# Representation as richness values: sp_pixel
rastercB2<- rasterbase
rastercB2[B2_pixel$Id]<- B2_pixel$sp_pixel
plot(rastercB2, add=T)
writeRaster(rastercB2, "rastercB2.tif", overwrite= T)


### Mapa por ecosistema
aa<- B2_pixel %>%dplyr:: group_by(Id) %>% dplyr::summarise(n= n_distinct(W_Ecosystm)) # pixels unicos por ecosistema
rasterEcosystem<- rasterbase
rasterEcosystem[B2_pixel$Id]<- B2_pixel$gridcode

match_id<- dplyr::select(B2_pixel, c("Id", "gridcode")) %>% distinct()

# raster a poligono
poly_ecosystems<- rasterToPolygons(rasterEcosystem) %>% st_as_sf() %>% dplyr::filter(!is.na(layer)) %>% dplyr::mutate(gridcode= layer) %>% 
  as.data.frame() %>% list(match_id) %>% join_all() %>% list(B2_pixel) %>% join_all()  %>%
  dplyr::select(c("Id", "gridcode", "chao_ecosystem", "W_Ecosystm", "name_clean", "geometry")) %>%  distinct() %>% st_as_sf()

st_write(poly_ecosystems, "pixels_ecosystems_spv5.shp")

write.csv(st_drop_geometry(poly_ecosystems), "pixels_ecosystems_spv5.csv")



