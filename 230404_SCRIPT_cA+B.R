# Chapter I: Conservation proposal for the Espeletiinae subtribe (Asteraceae) based on the Important Plant Areas (IPAs) methodology
# IPA methodology (Diazgranados & Castellanos-Castro 2018)

rm(list = ls())

## Installation Packages
packages_list<-list("data.table", "sp", "sf", "fasterize","raster", "rgdal", "stars", "foreign", "terra", "cli", "vegan", "tidyverse", "scales", "plyr", "fdth", "reshape", "reshape2","ggnewscale","raster", "rgeos", "ggplot2", "tmap", "pryr","base", "maptools", "ggmap", "rnaturalearth", "rnaturalearthdata", "rstudioapi", "magrittr", "pbapply", "rredlist", "grid", "gdalUtilities", "rCAT")
lapply(packages_list, library, character.only = TRUE)

# El poligono se cortó con capa de DEM de 90 m desde los 1000m
# Spatial Information Organization (Raster) 10km va a ser más rápido. raster es una grilla. El shp se rasteriza para: (i) shp tienden a ser vectores=pesado, (ii) dificiles de cuantificar por vertices. Mejor raster, es un archivo cuadrado, pixeles, grilla en si mism
setwd("C:/Users/Clari Mora/Dropbox/PC/Documents")
dirfile<- "C:/RESULTADOS IPAs_2022/Poligono_Area Estudio_1000 m_UA 10x10km/UA_1000_10x10.shp"
# cargar archivo (llamar el polígono) y transformarlo en raster con coordenadas planas (pasando de geográficas)
spatial_file<- st_read(dirfile) %>% st_transform(3395)

# una vez cargado, se resteriza el poligono (raster base de trabajo). 
# Se crea Id nuevo por cada pixel de 10 x 10 km2
# Sirve paa hacer mapa espacial. Es el raster base
rasterbase<- raster(extent(spatial_file),crs= paste0("+init=epsg:", 3395), res= 10000 )
tname3<- tempfile(fileext = ".tif")
t_file<- writeStart(rasterbase, filename = tname3,  overwrite=T); writeStop(t_file);
gdalUtilities::gdal_rasterize(dirfile, tname3, burn=1, at=T)

# generar Id de pixel, de 0 al maximo número. Creación de máscará de área de estudio raster
raster_file<- raster(tname3) %>% {terra::mask(setValues(., seq(ncell(.))), .)} %>% rast

# Copia de la máscara en 4326. Es mas fácil transformar la máscara que todos los datos
area_4326<- raster(raster_file) %>% projectRaster(crs = st_crs(4326)$proj4string, method = "ngb")
cell_area<- terra::cells(raster_file)

poly =rasterToPolygons(raster(raster_file)) %>% st_as_sf() %>% st_transform (4326) 
names(poly) [1]="Id"
st_write(poly, "poligono_pixel_all_ve3.shp")

## Load Data Occ Espeletiinae - cargar datos
glob <- read.csv("C:/RESULTADOS IPAs_2022/221103_OCC_SRTM_Copernicus_TNC_Useful_Coord.csv", sep=",", header = TRUE, row.names = NULL, na.strings = "")

# Cruce del poligno versus registros. Asignacion del Id del pixel {is.na = para que no haya NAs}
# Se reemplada , por . y que se vuelvan númericos (las coordenadas)
# se crea un Id coord, puede ayudar para el conteo. Son 2 Ids 
glob <- glob %>% 
  mutate(LATITUDE= as.numeric(gsub(",", ".", LATITUDE)), LONGITUD= as.numeric(gsub(",", ".", LONGITUD)), id_coord= group_indices_(., .dots=c("LATITUDE", "LONGITUD")) ) %>%
  dplyr::filter(!is.na(LATITUDE) | !is.na(LONGITUD))

# Cruce de registros espacializados
# se deben relacionar las columnas con coordenadas
glob_spatial<- glob %>% 
  st_as_sf(coords = c("LONGITUD","LATITUDE"), crs = 4326)

# Visualización puntos
plot(glob_spatial[, "geometry"])

# Cruzar mascara con registros con sistema de coordenadas. Cuando hay cruces, se asigna Id. Pero aquellos puntos que no estén con la mascará, eliminar
# Archivo espacial. es un archivo espacial, con ID de cell
# Convierte en un data frame para usar como tabla, y crea una celda, con el Id pero también funciona para filtrar por mascara
glob_spatial_mask<- as.data.frame(glob_spatial) %>% mutate(Id= raster::extract(area_4326, glob_spatial) ) %>%
  dplyr::filter(!is.na(Id))  %>%  st_as_sf() %>% st_transform(3395)

# Creación data.frame para manejo de información, en tabla
glob_data<- as.data.frame(glob_spatial_mask)

# Conteo de registros por especie (TODOS occ), por pixel diferenciado por coord
# n_occ_sp: # occ de especie por pixel
# n_occ_sp_pixel: # occ de sp por pixel. Conteo sp, por pixel y cantidad por pixel (suma de occ / pixel)
# se asume que cada registro es una coord diferente.:n_distinct: coord diferentes
data_id_sp<- glob_data %>% group_by(Id,name_clean) %>% dplyr::summarise(n_occ_sp_pixel = n_distinct(id_coord))

# Glob data agrupado por nombre de sp, cuantas veces esta la sp en total y por pixel unidos por name_clean
# sp por Id / pixel, y occ tales sp / Id. #  occ de la sp: cuantas sp hay en total
data_sp<- data_id_sp %>% group_by(name_clean) %>% dplyr::summarise(n_occ_sp_total = sum(n_occ_sp_pixel))

# Lista de registros totales tanto de ocurrencias, numero de especies por Id/pixel [base análisis]
data_summ<- list(data_id_sp, data_sp) %>% join_all()

# Exportar resultado
write.table(data_summ,file="occ_sp_pixel_v2.csv", sep = ",", row.names = TRUE, col.names=TRUE)

##############################################################################
## CRITERIO cA1: sp GLOBALLY THREATENED
# El área contiene poblaciones importantes de taxones amenazados globalmente
# Condicion: al menos el 1% de todos los especímenes registrados por especie amenazada a nivel global se encuentran dentro de un pixel
##############################################################################

# Filtrar por columna de amenazada, bajo columnas de IUCN. Crear vector único con valores requeridos
# dplyr::filter: que me filtre de la BD grande, las columnas que necesito (con %in% valores que quiero que coja) y creo un vector con las columnas que necesito
# da el valor de cuantas sp estan amenazadas
glob_sp_IUCN<- dplyr::filter(glob_data, IUCN_RL_CO %in% c("VU", "EN", "CR"))$name_clean %>% unique()

# Conteos bajo Thresholds por especie sobre el total de pixeles ocupados
# Cuántas occ hay por sp en toda la mascará y luego carlcular que % (1%) representa y si clasifican de acuerdo con el umbral
glob_data_IUCN_1perc<- dplyr::filter(data_summ, name_clean %in% glob_sp_IUCN) %>% 
  mutate(percent= n_occ_sp_pixel/ n_occ_sp_total) %>% ## add nueva columna con el calculo con el % calculado del umbral[mutate]
  dplyr::filter(percent>=0.01)

# por medio de un gropup_by, por pixel Id, cuántas sp con al menos el 1% cumplen para generar un map de distribucion
# Arroja cada sp
# n_distict: cuantas sp diferentes hay en el pixel. cA1 # de especpies por pixel.n_occ_sp_pixel: de la especie cuantas veces esta por pixel
glob_data_IUCN_1perc_spatial <- glob_data_IUCN_1perc %>% group_by(Id) %>% 
  dplyr::summarise(cA1= n_distinct(name_clean), name_clean=name_clean) %>% 
  distinct() %>% list(dplyr::select(data_summ, -n_occ_sp_total)) %>% join_all

# Creacion de tabla con columnas Id pixel, Criterio cA1 y occ por pixel.
glob_data_IUCN_1perc_spatial_uniques <- dplyr::select(glob_data_IUCN_1perc_spatial, -name_clean) %>% distinct()

# Mapeo de resultado y exportación
rastercA1<- rasterbase
rastercA1[glob_data_IUCN_1perc_spatial_uniques$Id]<- glob_data_IUCN_1perc_spatial_uniques$cA1
plot(raster_file, col= "gray")
plot(rastercA1, add=T)
names(rastercA1) = "cA1"

# organización información, por pixel y dentro de este cuántas sp y cuántas veces se encuentran. 
# Columna cA1_v1: cantidad de sp que caen en pixel. 
summ_sp_pixel_cA1 <- split(glob_data_IUCN_1perc_spatial, glob_data_IUCN_1perc_spatial$Id) %>% 
  lapply(function(x) {dplyr::mutate( x[1,], name_clean= paste(x$name_clean, collapse= ";" ) )}) %>%
  rbind.fill()

# Consolidar polígono (shp) y base en csv con información - resultado cA1
poly_cA1 =rasterToPolygons(rastercA1) %>% st_as_sf() %>% 
  as.data.frame() %>% dplyr::mutate(Id = terra::cells(rast(rastercA1))) %>%
  list(summ_sp_pixel_cA1) %>% join_all() %>% st_as_sf %>%
  st_transform(4326)

st_write(poly_cA1, "cA1_evaluado.shp")
write.table(glob_data_IUCN_1perc_spatial,file="cA1_evaluado.csv", sep = ",", row.names = TRUE, col.names=TRUE)

# extracción de datos con base en las columnas para análisis
# se filtra por la columna de conteos
# drop=TRUE: donde no hay nombres, se elimina
# Visializacion de columna ID pixel por sp hacia la derecha
draft_table_cA1 <- dcast(glob_data_IUCN_1perc_spatial, Id ~ name_clean, drop=TRUE, fill = 0, value.var = "n_occ_sp_pixel")
names(draft_table_cA1)
write.table(draft_table_cA1, file="ccA1_organice.csv", sep = ",", row.names = TRUE, col.names=TRUE)

# Organización y cálculo de Porcentajes de tabla organizada por Id y occ por especie, por cada sp en cada pixel
# Porcentaje, organizacion de los % (altos y bajos:Sort by %), normalizacion de datos (escala 0 a 1)
# Id: número de pixel en area de estudio
# name_clean: nombre sp
# perc_occ: porcentaje umbral (al menos 1%, es decir, igual o mayor)
# Cumple: si cumple con umbral o no
# Rank cA1: de cada especie, se saca su total: si cae en un 1 pixel pero esta en 50 se hace una relacion de %. Sobre este % se saca un ranking de los pixeles que más tienen al que menos (escala de 1 en adelante). De este rank, se hace una normalizacion
# a menor puntaje, más importancia de ese pixel.
# Norm cA1: corresponde a la organización de la columna del ranking, normalizada de llevada a escala de 0 a 1

percent_IUCN_1perc <- glob_data_IUCN_1perc %>% group_by(Id,name_clean) %>%
  dplyr::summarise(perc_occ = n_occ_sp_pixel/ n_occ_sp_total, threshold = n_occ_sp_total * 0.01) %>%
  mutate(cumple = ifelse (perc_occ >= 0.01, "Yes", "No")) %>%
  split(. $name_clean) %>% 
  lapply(function(x) mutate (x, RankcA1 = group_indices_(x, .dots=c("perc_occ"))) %>% mutate (RankcA1 = max(. $RankcA1) -RankcA1+1) %>% 
           mutate(NormcA1 = ((RankcA1 *100) / sum (. $RankcA1)) / 100)
           ) %>% rbind.fill()

write.table(percent_IUCN_1perc, file="cA1_1perc.csv", sep = ",", row.names = TRUE, col.names=TRUE)

##############################################################################
## CRITERION cA3: sp NATIONALLY THREATENED
# El área contiene poblaciones importantes de taxones amenazados nacionalmente
# Condicion: al menos el 10% de todos los especímenes registrados por especie amenazada en el país, se encuentran dentro de un pixel
##############################################################################

# Filtrar por columna de amenazada, por país por valor de IUCN. Crear vector único con valores requeridos
# a cada elemento por pais, se filtra por categoria

country_sp_IUCN <- split(glob_data, glob_data$PAIS) %>% 
  lapply(function(x) dplyr::filter(x, GEPC_Natio %in% c("VU", "EN", "CR"))$name_clean) %>%
  unlist() %>% unique()

# saber cuáles son las sp que de lo global, no están en lo nacional (es informativo)
sp_global_no_country <- glob_sp_IUCN %>% {.[!. %in% country_sp_IUCN]}
sp_coutry_no_global <- country_sp_IUCN %>% {.[!. %in% glob_sp_IUCN]}

# Conteos bajo Thresholds por especie sobre el total de pixeles ocupados
# Cuántas occ hay por sp en toda la mascara y luego carlcular que % (10%) representa y si clasifican de acuerdo con el umbral
country_data_IUCN_10perc <- dplyr::filter(data_summ, name_clean %in% country_sp_IUCN) %>% 
  mutate(percent= n_occ_sp_pixel/ n_occ_sp_total) %>%
  dplyr::filter(percent>=0.1)

# por medio de un gropup_by, por las columnas Id, cuántas sp con al menos el 1% cumplen para generar un map de distribucion
country_data_IUCN_10perc_spatial <- country_data_IUCN_10perc %>% group_by(Id) %>% 
  dplyr::summarise(cA3= n_distinct(name_clean), name_clean=name_clean) %>% 
  distinct() %>% list(dplyr::select(data_summ, - n_occ_sp_total)) %>% join_all

# Creacion de tabla con columnas Id pixel, Criterio cA3 y occ por pixel.
country_data_IUCN_10perc_spatial_uniques <- dplyr::select(country_data_IUCN_10perc_spatial, -name_clean) %>% distinct()

# Mapeo de resultado y exportación
rastercA3<- rasterbase
rastercA3[country_data_IUCN_10perc_spatial_uniques$Id]<- country_data_IUCN_10perc_spatial_uniques$cA3
plot (raster_file, col= "gray")
plot(rastercA3)
names(rastercA3) = "cA3"

# organización información, por pixel y dentro de este cuántas sp y cuántas veces se encuentran.
# Columna cA3: cantidad de sp que caen en pixel.
summ_sp_pixel_cA3 <- split(country_data_IUCN_10perc_spatial, country_data_IUCN_10perc_spatial$Id) %>% 
  lapply(function(x) {dplyr::mutate( x[1,], name_clean= paste(x$name_clean, collapse= ";" ) )})  %>%
  rbind.fill()

# Consolidar polígono (shp) y base en csv con información - resultado cA3
poly_cA3 =rasterToPolygons(rastercA3) %>% st_as_sf() %>% 
  as.data.frame() %>% dplyr::mutate(Id = terra::cells(rast(rastercA3))) %>%
  list(summ_sp_pixel_cA3) %>% join_all() %>% st_as_sf %>%
  st_transform(4326)

# Exportar
st_write(poly_cA3, "cA3_evaluado.shp")
write.table(country_data_IUCN_10perc_spatial,file="ccA3_evaluado.csv", sep = ",", row.names = TRUE, col.names=TRUE)

# extracción de datos con base en las columnas para análisis
draft_table_cA3 <- dcast(country_data_IUCN_10perc_spatial, Id ~ name_clean, drop=TRUE, fill = 0, value.var = "n_occ_sp_pixel")
write.table(draft_table_cA3, file="ccA3_organice.csv", sep = ",", row.names = TRUE, col.names=TRUE)

# Ogranización y cálculo de Porcentajes de tabla organizada por Id y occ por especie, por cada sp en cada pixel
# Porcentaje, organizacion de los % (altos y bajos:Sort by %), normalizacion de datos (escala 0 a 1)
# ultima columna con normalización de los datos
percent_country_IUCN_10perc <- country_data_IUCN_10perc %>% group_by(Id,name_clean) %>%
  dplyr::summarise(perc_occ_country = n_occ_sp_pixel/ n_occ_sp_total, threshold = n_occ_sp_total * 0.1) %>%
  mutate(cumple = ifelse (perc_occ_country >= 0.1, "Yes", "No")) %>%
  split(. $name_clean) %>% 
  lapply(function(x) mutate (x, RankcA3 = group_indices_(x, .dots=c("perc_occ_country"))) %>% mutate (RankcA3 = max(. $RankcA3) -RankcA3+1) %>% 
           mutate(NormcA3 = ((RankcA3 *100) / sum (. $RankcA3)) / 100)
  ) %>% rbind.fill()

write.table(percent_country_IUCN_10perc, file="cA3_10perc.csv", sep = ",", row.names = TRUE, col.names=TRUE)

##############################################################################
##### IUCN EVALUATION - EOO (extension de presencia) y AOO (área de ocupación)
#############################################################################

# BD insumo inical [glob <- read.csv(system.file("extdata", "221103_OCC_SRTM_Copernicus_TNC_Useful_Coord.csv")] filtrar por columna de sp y coord
occ <-data.frame(taxa=glob$name_clean, lat=glob$LATITUDE,long=glob$LONGITUD)

# Resultado de calculo de EOO y AOO, asignacion de categoria IUCN
resultsdfCat_IUCN <- ConBatch(occ,cellsize = 2000,"km", FALSE)

##############################################################################
##### CRITERIO cA4: occ EAR
# El área contiene poblaciones importantes de taxones endémicos (o casi endémicos) de rango altamente restringico (EAR)
# Condiciones: (i) Una especie EAR es aquella cuya extensión de ocurrencia es igual o menor a 100 km2; (ii) al menos el 10% de cada especie en rango EAR se encuentra en un pixel
##############################################################################

EOO_EAR = resultsdfCat_IUCN %>% dplyr::mutate(EOOkm2 = as.numeric(EOOkm2))
sp_EAR <- EOO_EAR %>% dplyr::filter (EOOkm2 <= 100)
View(sp_EAR)

# Calculo de cuántas hay por pixel y filtrarlo por columna de sp. y con la condicion de EAR
# SUMM= RESUMEN
sp_EAR_summ <- dplyr::filter(data_summ, name_clean %in% sp_EAR$taxa)
write.table(sp_EAR_summ,file="sp_EAR_cA4.csv", sep = ",", row.names = TRUE, col.names=TRUE)

# Calculo de umbral. Sacar frecuencias relativas al menor el 10% de todas las occ de sp
sp_EAR_10perc <- sp_EAR_summ %>% group_by(Id) %>% 
  mutate(percent = n_occ_sp_pixel / n_occ_sp_total) %>% 
  dplyr::filter(percent >= 0.1)
view(sp_EAR_10perc)

# Pixel, # occ, % que representan
# Se añade un filtro extra: cA4
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

# Ogranización y cálculo de Porcentajes de tabla organizada por Id y occ por especie, por cada sp en cada pixel
# Porcentaje, organizacion de los % (altos y bajos:Sort by %), normalizacion de datos (escala 0 a 1)

percent_cA4_EAR <- sp_EAR_10perc_pixel %>% group_by(Id,name_clean) %>%
  dplyr::summarise(perc_occ_EAR = n_occ_sp_pixel/ n_occ_sp_total, threshold = n_occ_sp_total * 0.1) %>%
  mutate(cumple = ifelse (perc_occ_EAR >= 0.1, "Yes", "No")) %>%
  split(. $name_clean) %>% 
  lapply(function(x) mutate (x, RankcA4 = group_indices_(x, .dots=c("perc_occ_EAR"))) %>% mutate (RankcA4 = max(. $RankcA4) -RankcA4+1) %>% 
           mutate(NormcA4 = ((RankcA4 *100) / sum (. $RankcA4)) / 100)
  ) %>% rbind.fill()

write.table(percent_cA4_EAR, file="cA4_EAR_result.csv", sep = ",", row.names = TRUE, col.names=TRUE)


##############################################################################
##### CRITERIO cA5: occ ERR
# El área contiene poblaciones importantes de taxones endémicos (o casi endémicos) de rango restringico (ERR)
# Condiciones: (i) Una especie ERR es aquella cuya extensión de ocurrencia es mayor de 100 km2 pero igual o menor de 5000 km2; (ii) al menos el 10% de cada especie en rango ERR se encuentra en un pixel
##############################################################################

EOO_ERR = resultsdfCat_IUCN %>% dplyr::mutate ( ( (EOOkm2 > 100) & (EOOkm2 <= 5000) ) )
sp_ERR <- EOO_EAR %>% dplyr::filter (EOOkm2 <= 100)
View(sp_ERR)

# Calculo de cuántas hay por pixel y filtrarlo por columna de sp. y con la condicion de EAR
# SUMM= RESUMEN
sp_ERR_summ <- dplyr::filter(data_summ, name_clean %in% sp_ERR$taxa)

# n_ERR = unique(sp_ERR_data_summ$name_clean) %>% length()
sp_ERR_occ = unique(sp_ERR$name_clean) %>% length()

# Calculo de umbral. Sacar frecuencias relativas al menor el 10% de todas las occ de sp
sp_ERR_10perc <- sp_ERR_summ %>% group_by(Id) %>% 
  mutate(percent= n_occ_sp_pixel/n_occ_sp_total) %>% 
  dplyr::filter(percent >= 0.1)
view(sp_ERR_10perc)

# Pixel, # occ, % que representan
# Se añade un filtro extra: cA5
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

# Ogranización y cálculo de Porcentajes de tabla organizada por Id y occ por especie, por cada sp en cada pixel
# Porcentaje, organizacion de los % (altos y bajos:Sort by %), normalizacion de datos (escala 0 a 1)

percent_cA5_ERR <- sp_ERR_10perc_pixel %>% group_by(Id,name_clean) %>%
  dplyr::summarise(perc_occ_ERR = n_occ_sp_pixel/ n_occ_sp_total, threshold = n_occ_sp_total * 0.1) %>%
  mutate(cumple = ifelse (perc_occ_ERR >= 0.1, "Yes", "No")) %>%
  split(. $name_clean) %>% 
  lapply(function(x) mutate (x, RankcA5 = group_indices_(x, .dots=c("perc_occ_ERR"))) %>% mutate (RankcA5 = max(. $RankcA5) -RankcA5+1) %>% 
           mutate(NormcA5 = ((RankcA5 *100) / sum (. $RankcA5)) / 100)
  ) %>% rbind.fill()

write.table(percent_cA5_ERR, file="cA5_ERR_result.csv", sep = ",", row.names = TRUE, col.names=TRUE)

########################
# CRITERIO B
# B1:El área tiene la más alta riqueza estimada de plantas y hongos nativos en ese tipo de ecosistema
#######################

#creación archivo temporal
tname_ecosystm <- tempfile(fileext = ".tif")

# se utiliza rasterbase creado previamente con los Id de pixel
t_file<- writeStart(rasterbase, filename = tname_ecosystm,  overwrite=T); writeStop(t_file);

# Importar capa de ecosistemas
# Informacion Ecosistema: Sayre, R., 2022, World Terrestrial Ecosystems (WTE) 2020
# U.S. Geological Survey data release, https://doi.org/10.5066/P9DO61LP.io
dir_ecosystems<- "C:/Users/Clari Mora/Dropbox/PC/Desktop/202210_BackUp_Info Clara/Doctorado/ESPELETIA/ArcGIS Insumos/230201_World Ecosystem SHP/World_Ecosystem.shp"

# rasterizar. No hay que convertir sistema de coord, rasterize lo hace
# gridcode: es el Id de capa de ecosistemas para unir con ecossitemas
gdalUtilities::gdal_rasterize(dir_ecosystems, tname_ecosystm,  at=T, a= "gridcode")

# reemplazar shp por dbd
library(foreign)
data_ecosystems<- read.dbf(gsub(".shp", ".dbf", dir_ecosystems))[, c("gridcode", "W_Ecosystm")] %>% distinct()

# rasterizacion de capa
rast_ecosystems<- rast(tname_ecosystm)

# convertir capa en tabla (disminuir gasto computacional)
# rast ecosystem: capa rasterizada
# terra:cell: que diga cuales celdas tienen valores - raster paquete base, terra segunda version de raster
# se debe poner gridcode: para que haga la unión 
# a que celda corresponde cada grid code (de mi sistema de pixeles)
# que llame las dos columnas de interés de mi sistema de grilla y ecosistemas. doy a pixeles el valor de gridcode. Id nunca cambia
table_ecosystems <- as.data.frame(rast_ecosystems) %>% setNames("gridcode") %>% mutate(Id= terra::cells(rast_ecosystems)) %>%
  list(data_ecosystems, data_summ) %>% join_all()

# quitar los NA de aquellas sp que que no caen en algun ecosistema
# sp por ecosistema
# nsp_ numero de especues por pixel
sp_ecosystems <-table_ecosystems %>% dplyr::filter(!is.na(name_clean)) %>% 
  group_by(W_Ecosystm) %>% dplyr::summarise(nsp_ecosystem= n_distinct(name_clean))

# calculo de cuantos occ hay por ecosistemas: cuantos distintas sp hay
sp_pixel <- table_ecosystems %>% group_by(Id) %>% dplyr::summarise(nsp_pixel= n_distinct(name_clean), W_Ecosystm=W_Ecosystm)

# listado ecosisema, # sp por pixel y por ecosistema: este es para hacer el join entr informacion de sp y ecosistemas
# # de sp por pixel y # sp por ecosistema
# tenemos riqueza real y proyectada y todos aquellos ecosistemas que tengan mas de 50% que se repoyectan : se aceptan (quitando 3, 2 1, registros)
# no se filta por un numero especifico, por que sería arbitrario
sp_pixel_ecosystem <- list(sp_ecosystems, sp_pixel) %>% join_all() %>% mutate(perc= nsp_pixel/ nsp_ecosystem) %>%
  dplyr::filter()

write.table(sp_pixel_ecosystem,file="sp_pixel_W_Ecosystm_v1.csv", sep = ",", row.names = TRUE, col.names=TRUE)

# que parta en listados de especices por ecosistema
# de los que tienen sp
# cuando se cargan cosas con factor, existen como 0, enton rescribir factor con los valores que si quedaron
#unique indica cuales son los ecosistemas
list_ecosystems <- dplyr::filter(table_ecosystems, W_Ecosystm %in% unique(sp_ecosystems$W_Ecosystm) ) %>% 
  dplyr::mutate(W_Ecosystm= factor(W_Ecosystm)) %>% 
  split(.$W_Ecosystm)

# riqueza estimada por ecosistema
# chao: cuantos podrian existir
# se borran los NAs
# SPECIES: REPRSENTATIVIDAD: DEBERIA DECIR CUANTO DA ESPECIES EN CHAO
# > 1 sobre especies, filtrando por que haya mas de 1 especie. Por que cuando chao daba igual a 1, daba NA en tabla
# chao siempre dará por ecosistema no por pixel
chao_ecosystems<- lapply(list_ecosystems, function(x){ print(unique(x$W_Ecosystm))
  
  tryCatch({
    
    draft_ecosystems <- dcast(x, Id ~ name_clean, drop=TRUE, fill = 0, value.var = "n_occ_sp_pixel")  %>%
      {.[, !names(.) %in% "NA"]} %>% column_to_rownames("Id")
    
    rich_ecosystems_estimated <- specpool(draft_ecosystems) %>% mutate(W_Ecosystm= unique(x$W_Ecosystm))
  }, error= function(e) {data.frame(W_Ecosystm=unique(x$W_Ecosystm))})
  
}) %>% rbind.fill() %>% mutate(repres= Species/chao) %>% dplyr::filter(Species>1)

write.csv(chao_ecosystems, file="sp_chao_ecosystems_v1.csv", sep = ",", row.names = TRUE, col.names=TRUE)   

###################
## CRITERIO B2
# B2:Aquella/s  área/s  con  riqueza  estimada  de  plantas  complementaria  a  B1,  de  manera  que  en  conjunto sumen el 10% de la biodiversidad total del ecosistema
# Condición cB2: cuando  la AIP seleccionada por el criterio B1 contenga el 10% o  más  de  la riqueza total estimada para el tipo de ecosistema, no será necesaria la aplicación del criterio B2.
##################

# columna repres: que representatividad tiene el pixel respecto al total: # sp / en todo el ecosistema = indice chao
# ceiling: para redondear decimales, sino da unidad completa
# todo lo que tenga mas del 10% será b2
# los que tienen NAs no los está calculando. Se borra de name_clean no se tendrá en cuenta
chao_richness<- data.frame(chao_ecosystem= chao_ecosystems$chao, W_Ecosystm= chao_ecosystems$W_Ecosystm) %>%
  list(table_ecosystems) %>% join_all() %>% group_by(Id) %>%
  dplyr::filter(!is.na(name_clean)) %>% 
  dplyr::mutate(sp_pixel= n_distinct(name_clean)) %>% dplyr::rowwise() %>%
  dplyr::mutate(chao_ecosystem= ceiling(chao_ecosystem)) %>% dplyr::mutate(repres=sp_pixel / (chao_ecosystem)) %>% 
  dplyr::mutate(B2= ifelse(repres>=0.1, "SI", NA))


# pixels con representatividad > del 10%
# visualiza tabla con chao, Id, name_clean, occ, sp_pixel, representatividad de acuerdo con chao, y columna de si cumple B2
B2_pixel<- dplyr::filter(chao_richness, !is.na(B2)) %>%  data.frame()
write.csv(B2_pixel, file = "cB2_complementario.csv")

# Salida de mapa cB2, los complementarios
# parto de una riqueza estimada por ecosistema, se indica que esos son los pixeles más importantes por ecosistma por los datos que tienen
# Representación como valores de riqueza: sp_pixel
rastercB2<- rasterbase
rastercB2[B2_pixel$Id]<- B2_pixel$sp_pixel
plot(rastercB2)
writeRaster(rastercB2, "rastercB2.tif")

# para sacar los mapas
# stars: potearlo mejor. Para mejorar visualización
stars_B2 <- stars::st_as_stars(rastercB2)
data_B2<- as.data.frame(stars_B2) %>% dplyr::filter(!is.na(layer)) 

# Exportar como puntos 
points_B2<- data_B2 %>% mutate(lat= y , long= x) %>%  st_as_sf(coords =c("long", "lat"), crs = 3395) %>% st_transform(4326)
st_write(points_B2, "points_B2.shp")


# Mapa de salida con resultados de cB2
  gg_Figb2<- ggplot() + geom_stars(data = stars_altitude, na.rm=T) + 
  scale_fill_gradient2(low = "white", high = "gray10", mid="gray", midpoint =3000, na.value = NA, guide = "none") + 
  geom_polygon(data = join_adm2, aes(x = long, y = lat, group = group), color = "gray0", fill = NA) + 
  new_scale_fill() +
  geom_point(data= data_B2, aes(x=x, y=y,color=layer, size=layer), alpha= 0.5, shape= 15)+
  scale_color_distiller(palette= "Spectral",breaks = seq(max(data_B2$layer)), "Número de\ndetonantes",  guide = guide_legend())+
  scale_size_continuous(breaks = seq(max(data_B2$layer)), "Número de\ndetonantes")+
  coord_fixed(xlim = ext_Andes[1:2], ylim = ext_Andes[3:4] )+
  theme_void()+
  theme(legend.position = "bottom", text = element_text(size = 10))+
  labs(x = "Longitude", y = "Latitude", title = "Número de detonantes\ncB2 por IPA potencial")

write.table(chao_richness, file="chao_richness_W_Ecosystm_v1.csv", sep = ",", row.names = TRUE, col.names=TRUE)
write.csv(sp_ecosystems, file="sp_ecosystems_W_Ecosystm_v1.csv", sep = ",", row.names = TRUE, col.names=TRUE)                                                                    

