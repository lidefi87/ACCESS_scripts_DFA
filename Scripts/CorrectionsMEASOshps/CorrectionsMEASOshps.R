############################################################################
#Obtaining MEASO limits from measoshapes package
#https://australianantarcticdivision.github.io/measoshapes/
#
#Shapefiles are currently available within the package, but need some work
#before they can be exported and used outside R.
############################################################################

# Loading libraries -------------------------------------------------------
library(measoshapes)
library(sf)
library(tidyverse)


# Accessing names of MEASO regions ----------------------------------------
#Names of measo regions
measo_reg_names <- measo_names %>% 
  #Include the zone name for all areas with names ending in T
  mutate(zone = case_when(str_detect(name, "T$") ~ "Temperate",
                          TRUE ~ zone))


# Adding supporting info to shapefiles ------------------------------------
#South Polar Stereographic
#Adding supporting information to MEASO regions shapefile
measo <- measo_regions05 %>% 
  #Only using unique values to avoid duplication
  left_join(measo_reg_names %>% distinct(), by = "name")

#Check crs
st_crs(measo)

#Plot MEASO regions
ggplot(measo, aes(fill = zone))+
  geom_sf(col = NA)

#Save shapefile to disk - append parameter set to F so it overwrites shapefile
st_write(measo, "Outputs/measo.shp", append = F)


#WGA84
#Adding supporting information to MEASO regions shapefile
measo_wgs84 <- measo_regions05_ll %>% 
  #Only using unique values to avoid duplication
  left_join(measo_reg_names %>% distinct(), by = "name")

#Check crs
st_crs(measo_wgs84)

#Plot MEASO regions
ggplot(measo_wgs84, aes(fill = zone))+
  geom_sf()

#Save shapefile to disk - append parameter set to F so it overwrites shapefile
st_write(measo_wgs84, "Outputs/measo_wgs84.shp", append = F)


# Plot Southern Ocean -----------------------------------------------------
#Access coastline of Antarctica
library(ggOceanMapsData)
library(ggOceanMaps)

Antarctica <- st_read("~/GIS DataBase/CoastlineAntarctica/add_coastline_medium_res_polygon_v7_4.shp")
Antarctica <- Antarctica %>% 
  mutate(colour = case_when(str_detect(surface, ".*ice.*") ~ "#ffffff",
                            TRUE ~ "#dfdfdf"))

world <- st_read("~/GIS DataBase/WorldCountries/Countries_30S_wgs84.shp")
world <- st_transform(world, st_crs(Antarctica))

world <- map_data("world") %>% 
  st_transform(world, st_crs(Antarctica))


measo <- st_transform(measo, st_crs(Antarctica))

#Make graph
basemap(limits = -30, land.col = "#dfdfdf", land.border.col = "#6f787f", bathymetry = T)+
  geom_sf(data = measo, aes(fill = fill))+
  geom_sf(data = Antarctica, aes(fill = colour), col = "#6f787f")+
  scale_fill_identity()

ggplot(measo, aes(fill = fill))+
  geom_sf(col = NA)+
  geom_sf(data = Antarctica, aes(fill = colour))+
  geom_sf(data = world, aes(fill = "#dfdfdf"))
  # scale_fill_identity()
  # theme_bw()+





SOmap(trim = -30)
SOplot(st_geometry(measo), border = NA)

