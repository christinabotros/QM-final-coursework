# Library packages which will be used to run the SIM and analysis 
library(sp)
library(reshape2)
library(geojsonio)
library(rgdal)
library(downloader)
library(maptools)
library(dplyr)
library(broom) 
library(stplanr)
library(ggplot2)
library(MASS)
library(sf)
library(tmap)
library(tmaptools)
library(stringr)
library(ggplot2)
library(leaflet)
library(tidyverse)
library(here)
library(downloader)
library(readxl)
library(janitor)
library(osrm)

# Loading in the LSOA shapefiles for Leeds. 

boundary <- "https://borders.ukdataservice.ac.uk/ukborders/servlet/download/dynamic/48AD826F4A472DADD9161099874237012/16109987426482818050824586534527/BoundaryData.zip" 
#download from the URL
download(boundary, dest="dataset.zip", mode="wb")
#unzip into a new directory called data
unzip("dataset.zip",exdir="./Data")
#get the file names from within the zip file
filename <- list.files("./Data")
#in order to read the OA boundary please clone/download the BoundaryData folder from the github repro - everything else is accessible via https links
here::here()
LSOA <- st_read(here::here("Data", "england_lsoa_2011.shp")) %>%
  st_transform(27700)

LSOA %>% 
  plot() 

## Geolytix Retail Points dataset 

Geolytix update their Retail Points  dataset several times a year and this research uses the latest version (August 2020). Within Geolytix grocery retailers’ floorspace is expressed as a size band (e.g., ‘between 3,013 and 15,069 sq.ft.’) rather than an exact number, which is why GMAP (2016) was used collaboratively to obtain the most precise and accurate profiles.  

# previous years grocery points available from here <- "https://drive.google.com/file/d/1B8M7m86rQg2sx2TsHhFa2d-x-dZ1DbSy/view?usp=sharing"

retail_points <- read.csv("https://raw.githubusercontent.com/christinabotros/Leeds-Food-Desert-Analysis/main/geolytix_retailpoints_v17_202008.csv")   
head(retail_points)
plot(retail_points$bng_e, retail_points$bng_n)

retail_sf <- st_as_sf(retail_points, coords = c("bng_e", "bng_n"), 
                      crs = 27700)

# Clip the retail points to only be for those in Leeds
retail_leeds <- retail_sf[LSOA,]


retailers_map <- tm_shape(LSOA$geometry) + 
  tm_polygons(col="white", alpha = 0.3) + 
  tm_shape(retail_leeds) + 
  tm_symbols(col ="red", scale =.3)
retailers_map

retail_leeds <- retail_leeds %>% rename('Grocery Retailers' = retailer)

finalmap <- tm_shape(LSOA) + tm_borders(col="black") +
  tm_shape(retail_leeds) + tm_symbols(col="Grocery Retailers", scale = 0.5, palette = "Accent") 
finalmap

tmap_save(finalmap, "LSOA grocery retailers.png")

#write.csv(retail_leeds, "retail_leeds.csv")

# unique names of the retailers

retailers <- unique(retail_leeds$retailer.x)

# [1] "Aldi"                   "Asda"                  
# [3] "Marks and Spencer"      "Morrisons"             
# [5] "Sainsburys"             "Tesco"                 
# [7] "Waitrose"               "Lidl"                  
# [9] "Costco"                 "Makro"                 
# [11] "The Co-operative Group" "Iceland"               
# [13] "Farmfoods"              "Heron"                 
# [15] "Jack Fultons"           "Spar

```

# 6.Distance Matrix 

# Calculating the distance between every OA and grocery store point. We have 2453 OAs and 205 retailers. 

retail_leeds_pts <- st_centroid(retail_leeds)
retail_leeds_pts <- retail_leeds_pts[,1]
qtm(retail_leeds_pts)
lsoa_leeds_pts <- st_centroid(LSOA)
lsoa_leeds_pts <- lsoa_leeds_pts[,3] %>% rename(., id = code)
#qtm(oa_leeds_pts)


### Calculate distance matrix [Cij

#create a vector of all origin and destination points
all_od_pts <- rbind(lsoa_leeds_pts, retail_leeds_pts)
plot(all_od_pts)
summary(all_od_pts)

#create some vectors of the IDs
LSOA_Origin_codes <- LSOA$code #vector of all OA code names for your origins (that you can link back to the point geometries)
destination_shops <- retail_leeds$id #vector of codes for your shops (that you can link back to your shop point data)

#this should create a square(ish) matrix from the list of origin and destination codes
tb <- as_tibble(matrix(nrow = length(LSOA_Origin_codes), ncol = length(destination_shops), dimnames = (list(LSOA_Origin_codes,destination_shops))))

#OK, this is a proper mess, but it works - can be tidied afterwards. 
#The idea is to get the rownames to the left of the matrix. Can definitely 
#do with in one step with relocate in dplyr, but here we go...
#create a new column to store the row (origin) names
tb <- tb %>% 
  mutate(row_name = LSOA_Origin_codes) %>% 
  #and then turn this column into some actual 'row names' i.e. not a real column but names for the rows
  column_to_rownames(var = "row_name")

#now because I couldn't do it in one step, now create a new column from rownames
tb_1 <- tb %>% rownames_to_column(var = "orig")

#now pivot this longer into a new paired-list of origins and destinations
tb_long <- pivot_longer(tb_1, cols = 2:ncol(tb_1), names_to = "dest")
tb_long$value <- 1

#now generate some staight-line flow lines. We could try and route these along roads
#but given how many, this would totally break your computer. Start easy. 
travel_lines <- od2line(flow = tb_long, zones = all_od_pts, origin_code = "orig", dest_code = "dest")
travel_lines

#test to see if it's worked - don't try and plot the whole thing or R will cry!
#sub <- travel_lines[1:10000,]
#tmap_mode("view")
#qtm(sub)

#now calculate some distances
distance_matrix <- geo_length(travel_lines)
#now attach this back to travel_lines
travel_lines$dist <- distance_matrix

# 8.  Adding in actual floorspace from GMAP (2016)

## SUPPLY SIDE ### FLOORSPACE Wj
#attractiveness of a store is given by its size for this model

floorspace <- read.csv("https://raw.githubusercontent.com/christinabotros/Leeds-Food-Desert-Analysis/main/floorspace%20dataset.csv")

retail_leeds <- merge(retail_leeds, 
                      floorspace, 
                      by.x = "id", 
                      by.y = "id")


# 9. SPATIAL INTERACTION MODEL - Production-Constrained

# (1) $$T_{ij} = A_i O_i W_j^\alpha d_{ij}^-\beta $$
  
#where

#(2) $$O_{i} = \sum_j T_{ij}$$
  
#and 

#(3) $$ A_i = \frac{1}{\sum_j W_j^\alpha d_{ij}^-\beta}$$
  

# Adding in variables to the distance matrix
travel_lines1 <- travel_lines
travel_lines1$floorspace <- retail_leeds$Floorspace[match(travel_lines1$dest, retail_leeds$id)]
travel_lines1$retailers <- retail_leeds$retailer.x[match(travel_lines1$dest, retail_leeds$id)]

## Getting Ai 
# note we can adjust brand attractive so that market share within Leeds matches that within the 
#Kantar WorldPanel grocery market shares
alpha = 1  

#beta can be calibrated to reflect people's ability to travel different distances 
beta = -0.3

#calculate some new wj^alpha and dij^beta values
wj2_alpha <- travel_lines1$floorspace^alpha
dist_beta <- travel_lines1$dist^beta
#calculate the first stage of the Ai values
travel_lines1$Ai1 <- wj2_alpha*dist_beta

#HANSEN ACCESSIBILITY SCORE

hasen <- dplyr::select(travel_lines1, c(orig, dest, Ai1)) %>% st_drop_geometry()
hasenmatrix <- hasen %>% pivot_wider(names_from = dest, values_from =Ai1) 

write.csv(hasenmatrix, "hasenmatrix.csv")


#now do the sum over all js bit
Hansenn <- travel_lines1 %>% group_by(orig) %>% summarise(hansenscore = sum(Ai1))  %>% st_drop_geometry()

hansen_map <- LSOA %>%
  merge(., 
        Hansenn, 
        by.x = "code", 
        by.y = "orig")

# qtm(hansen_map, 
#     fill = "Hansen Accessibility Score")

hansenmapp <-  tm_shape(hansen_map) + tm_fill("hansenscore") + 
  tm_borders(col="light grey", alpha=0.2) + 
  tm_layout(title = "Hansen Accessibility Score in Leeds LSOAs", title.size = 1.1, title.position = c("centre", "top"))
hansenmapp 

tmap_save(hansenmapp, "LSOA hansen accessibility.png")

#Reading in mid-year population estimates from ONS.
url_lsoa_pop_zip <- "https://www.ons.gov.uk/file?uri=%2fpeoplepopulationandcommunity%2fpopulationandmigration%2fpopulationestimates%2fdatasets%2flowersuperoutputareamidyearpopulationestimates%2fmid2019sape22dt2/sape22dt2mid2019lsoasyoaestimatesunformatted.zip"
#download from the URL
download(url_lsoa_pop_zip, dest="dataset.zip", mode="wb")
#unzip into a new directory called data
unzip("dataset.zip",exdir="./Data")
#get the file names from within the zip file
filename <- list.files("./Data")
#read the sheet you want from the file
lsoa_pop <- read_xlsx(here::here("Data", "SAPE22DT10c-mid-2019-coa-unformatted-syoa-estimates-yorkshire-and-the-humber.xlsx"), sheet="Mid-2019 Persons", skip = 3, col_names = T)
lsoa_pop <- clean_names(oa_pop)

LSOA <- LSOA %>%
  merge(., 
        lsoa_pop, 
        by.x = "code", 
        by.y = "lsoa11cd")

Hansenn$total_pop <- LSOA$total_pop[match(Hansenn$orig, LSOA$code)]
Hansenn$provision <- Hansenn$hansenscore / Hansenn$total_pop

hansen_map <- LSOA %>%
  merge(., 
        Hansenn, 
        by.x = "code", 
        by.y = "orig")

hansenmap <-  tm_shape(hansen_map) + tm_fill("provision") + 
  tm_borders(col="light grey", alpha=0.2) + 
  tm_layout(title = "Hansen Accessibility Score in Leeds LSOAs", title.size = 1.1, title.position = c("centre", "top"))
hansenmap

write.csv(Hansenn, "hasen provision final.csv")
