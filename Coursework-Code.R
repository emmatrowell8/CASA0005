#load libraries
library(tidyverse)
library(sf)
library(dplyr)
library(tmap)
library(tmaptools)
library(spdep)
library(RColorBrewer)

#load data
nyc_community_districts <- st_read('Data/nycd_20d/nycd.shp')
nyc_pm2.5 <- read.csv('Data/Neighbourhood_means/PM2.5_2019.csv')
nyc_no2 <- read.csv('Data/Neighbourhood_means/NO2_2019.csv')
nyc_so2 <- read.csv('Data/Neighbourhood_means/SO2_2019.csv')
nyc_poverty <- read.csv('Data/poverty-data-2018.csv') 
nyc_race <- read.csv('Data/ethnicity-2010.csv')
nyc_income <- read.csv('Data/CCC Data Download_20210105_143328175/Median Incomes.csv')

#tidy pollution and income data
nyc_no2_tidy <- nyc_no2 %>%
  subset(select = -c(1,3,6:8,10:12)) 
nyc_no2_tidy<- nyc_no2_tidy %>%
  filter(GeoTypeName == 'Neighborhood (Community District)')

nyc_pm2.5_tidy <- nyc_pm2.5 %>%
  subset(select = -c(1,3,6:8,10:12)) 
nyc_pm2.5_tidy<- nyc_pm2.5_tidy %>%
  filter(GeoTypeName == 'Neighborhood (Community District)')

nyc_so2_tidy <- nyc_so2 %>%
  subset(select = -c(1,3,6:8,10:12)) 
nyc_so2_tidy<- nyc_so2_tidy %>%
  filter(GeoTypeName == 'Neighborhood (Community District)')

nyc_income_tidy <- nyc_income %>%
  filter(TimeFrame == '2019') %>%
  filter(Household.Type == 'All Households')
nyc_income_tidy <- nyc_income_tidy[-c(60:65),]

#summary statistics
summary(nyc_no2_tidy$Mean..ppb.)
summary(nyc_pm2.5_tidy$Mean..mcg.per.cubic.meter.)
summary(nyc_so2_tidy$Mean..ppb.)
summary(nyc_poverty$Percent_in_poverty)
summary(nyc_race$X..of.population.non.white)
summary(nyc_income_tidy$Median_Income)

# merge district code to race and poverty data 
nyc_CD_IDs <- read.csv('Data/CD-Geography-IDs.csv') 
nyc_CD_IDs <- nyc_CD_IDs[c(1,2)]
nyc_CD_IDs <- nyc_CD_IDs[-c(60:62),]

nyc_community_districts <- merge(nyc_community_districts, nyc_CD_IDs, by.x = 'BoroCD', by.y = 'Geography.ID')

#join data to shapefile 
CD_no2 <- merge(nyc_community_districts, nyc_no2_tidy, by.x = 'BoroCD', by.y = 'Geography.ID')
CD_pm2.5 <- merge(nyc_community_districts, nyc_pm2.5_tidy, by.x = 'BoroCD', by.y = 'Geography.ID')
CD_so2 <- merge(nyc_community_districts, nyc_so2_tidy, by.x = 'BoroCD', by.y = 'Geography.ID')
CD_poverty <- merge(nyc_community_districts, nyc_poverty, by.x = 'BoroCD', by.y = 'Geography.ID')
CD_race <- merge(nyc_community_districts, nyc_race, by.x = 'BoroCD', by.y = 'Geography.ID')
CD_income <- merge(nyc_community_districts, nyc_income_tidy, by.x = 'BoroCD', by.y = 'Fips')

#st as sf
CD_no2_sf <- st_as_sf(CD_no2)
CD_pm2.5_sf <- st_as_sf(CD_pm2.5)
CD_so2_sf <- st_as_sf(CD_so2)
CD_poverty_sf <- st_as_sf(CD_poverty)
CD_race_sf <- st_as_sf(CD_race)
CD_income_sf <- st_as_sf(CD_income)

#analysis - spatial autocorrelation
#neighbours
nyc_neighbours <- poly2nb(nyc_community_districts)
plot(nyc_community_districts$geometry)
plot(nyc_neighbours, coords = nyc_community_districts$geometry, add  = TRUE, col = 'red')

#spatial weights matrix
spatial_weights <- nb2listw(nyc_neighbours, zero.policy = TRUE)

#Global moran's I
no2_global_moran <- moran.test(CD_no2_sf$Mean..ppb., spatial_weights, zero.policy = TRUE)
pm2.5_global_moran <- moran.test(CD_pm2.5_sf$Mean..mcg.per.cubic.meter., spatial_weights, zero.policy = TRUE)
so2_global_moran <- moran.test(CD_so2_sf$Mean..ppb., spatial_weights, zero.policy = TRUE)
poverty_global_moran <- moran.test(CD_poverty_sf$Percent_in_poverty, spatial_weights, zero.policy = TRUE)
race_global_moran <- moran.test(CD_race_sf$X..of.population.non.white, spatial_weights, zero.policy = TRUE)
income_global_moran <- moran.test(CD_income_sf$Median_Income, spatial_weights, zero.policy = TRUE)

#Results
no2_global_moran
pm2.5_global_moran
so2_global_moran
poverty_global_moran
race_global_moran
income_global_moran

#Local moran's I
no2_local_moran <- localmoran(CD_no2_sf$Mean..ppb., spatial_weights, zero.policy = TRUE)
pm2.5_local_moran <- localmoran(CD_pm2.5_sf$Mean..mcg.per.cubic.meter., spatial_weights, zero.policy = TRUE)
so2_local_moran <- localmoran(CD_so2_sf$Mean..ppb., spatial_weights, zero.policy = TRUE)
poverty_local_moran <- localmoran(CD_poverty_sf$Percent_in_poverty, spatial_weights, zero.policy = TRUE)
race_local_moran <- localmoran(CD_race_sf$X..of.population.non.white, spatial_weights, zero.policy = TRUE)
income_local_moran <- localmoran(CD_income_sf$Median_Income, spatial_weights, zero.policy = TRUE)

#Results
no2_local_moran
pm2.5_local_moran
so2_local_moran
poverty_local_moran
race_local_moran
income_local_moran

#plot local moran's I
no2_local_moran_map <- cbind(nyc_community_districts, no2_local_moran)
names(no2_local_moran_map)[5] <- 'Local Morans I'

pm2.5_local_moran_map <- cbind(nyc_community_districts, pm2.5_local_moran)
names(pm2.5_local_moran_map)[5] <- 'Local Morans I'

so2_local_moran_map <- cbind(nyc_community_districts, so2_local_moran)
names(so2_local_moran_map)[5] <- 'Local Morans I'

poverty_local_moran_map <- cbind(nyc_community_districts, poverty_local_moran)
names(poverty_local_moran_map)[5] <- 'Local Morans I'

race_local_moran_map <- cbind(nyc_community_districts, race_local_moran)
names(race_local_moran_map)[5] <- 'Local Morans I'

income_local_moran_map <- cbind(nyc_community_districts, income_local_moran)
names(income_local_moran_map)[5] <- 'Local Morans I'

#plots of local moran's I 
tm_shape(no2_local_moran_map) +
  tm_fill(col = "Local Morans I",
          style = "pretty",
          palette = 'PRGn') +
  tm_borders(alpha = 0.2) +
  tm_layout(main.title = 'Local Morans I of NO2 concentration',
            main.title.position = 'centre',
            frame = FALSE) +
  tm_compass(position = c('right', 'bottom')) +
  tm_scale_bar(position = c('right', 'bottom'), width = 0.3)

tm_shape(pm2.5_local_moran_map) +
  tm_fill(col = "Local Morans I",
          style = "pretty",
          palette = 'PRGn') +
  tm_borders(alpha = 0.2) +
  tm_layout(main.title = 'Local Morans I of PM2.5 concentration',
            main.title.position = 'centre',
            frame = FALSE) +
  tm_compass(position = c('right', 'bottom')) +
  tm_scale_bar(position = c('right', 'bottom'), width = 0.3)

tm_shape(so2_local_moran_map) +
  tm_fill(col = "Local Morans I",
          style = "pretty",
          palette = 'PRGn') +
  tm_borders(alpha = 0.2) +
  tm_layout(main.title = 'Local Morans I of SO2 concentration',
            main.title.position = 'centre',
            frame = FALSE) +
  tm_compass(position = c('right', 'bottom')) +
  tm_scale_bar(position = c('right', 'bottom'), width = 0.3)

tm_shape(poverty_local_moran_map) +
  tm_fill(col = "Local Morans I",
          style = "pretty",
          palette = 'PRGn') +
  tm_borders(alpha = 0.2) +
  tm_layout(main.title = 'Local Morans I of poverty distribution',
            main.title.position = 'centre',
            frame = FALSE) +
  tm_compass(position = c('right', 'bottom')) +
  tm_scale_bar(position = c('right', 'bottom'), width = 0.3)

tm_shape(race_local_moran_map) +
  tm_fill(col = "Local Morans I",
          style = "pretty",
          palette = 'PRGn') +
  tm_borders(alpha = 0.2) +
  tm_layout(main.title = 'Local Morans I of racial distribution',
            main.title.position = 'centre',
            frame = FALSE) +
  tm_compass(position = c('right', 'bottom')) +
  tm_scale_bar(position = c('right', 'bottom'), width = 0.3)

tm_shape(income_local_moran_map) +
  tm_fill(col = "Local Morans I",
          style = "pretty",
          palette = 'PRGn') +
  tm_borders(alpha = 0.2) +
  tm_layout(main.title = 'Local Morans I of income distribution',
            main.title.position = 'centre',
            frame = FALSE) +
  tm_compass(position = c('right', 'bottom')) +
  tm_scale_bar(position = c('right', 'bottom'), width = 0.3)

#Getis Ord Gi* 
Gi_no2 <- localG(CD_no2_sf$Mean..ppb., spatial_weights, zero.policy = TRUE)
CD_Gi_no2 <- cbind(nyc_community_districts, as.matrix(Gi_no2))
names(CD_Gi_no2)[5] <- 'Getis Ord Gi*'

tm_shape(CD_Gi_no2) +
  tm_fill("Getis Ord Gi*",
          palette = '-RdBu',
          style = 'pretty') +
  tm_borders(alpha = 0.2) +
  tm_layout(main.title = "Hot and cold spots of NO2 Concentration",
            main.title.position = 'centre',
            frame = FALSE) +
  tm_compass(position = c('right', 'bottom')) +
  tm_scale_bar(position = c('right', 'bottom'), width = 0.3)

Gi_pm2.5 <- localG(CD_pm2.5_sf$Mean..mcg.per.cubic.meter., spatial_weights, zero.policy = TRUE)
CD_Gi_pm2.5 <- cbind(nyc_community_districts, as.matrix(Gi_pm2.5))
names(CD_Gi_pm2.5)[5] <- 'Getis Ord Gi*'

tm_shape(CD_Gi_pm2.5) +
  tm_fill("Getis Ord Gi*",
          palette = '-RdBu',
          style = 'pretty') +
  tm_borders(alpha = 0.2) +
  tm_layout(main.title = "Hot and cold spots of PM2.5 Concentration",
            main.title.position = 'centre',
            frame = FALSE) +
  tm_compass(position = c('right', 'bottom')) +
  tm_scale_bar(position = c('right', 'bottom'), width = 0.3)

Gi_so2 <- localG(CD_so2_sf$Mean..ppb., spatial_weights, zero.policy = TRUE)
CD_Gi_so2 <- cbind(nyc_community_districts, as.matrix(Gi_so2))
names(CD_Gi_so2)[5] <- 'Getis Ord Gi*'

tm_shape(CD_Gi_so2) +
  tm_fill("Getis Ord Gi*",
          palette = '-RdBu',
          style = 'pretty') +
  tm_borders(alpha = 0.2) +
  tm_layout(main.title = "Hot and cold spots of SO2 Concentration",
            main.title.position = 'centre',
            frame = FALSE) +
  tm_compass(position = c('right', 'bottom')) +
  tm_scale_bar(position = c('right', 'bottom'), width = 0.3)

Gi_poverty <- localG(CD_poverty_sf$Percent_in_poverty, spatial_weights, zero.policy = TRUE)
CD_Gi_poverty <- cbind(nyc_community_districts, as.matrix(Gi_poverty))
names(CD_Gi_poverty)[5] <- 'Getis Ord Gi*'

tm_shape(CD_Gi_poverty) +
  tm_fill("Getis Ord Gi*",
          palette = '-RdBu',
          style = 'pretty') +
  tm_borders(alpha = 0.2) +
  tm_layout(main.title = "Hot and cold spots of Poverty distribution",
            main.title.position = 'centre',
            frame = FALSE) +
  tm_compass(position = c('right', 'bottom')) +
  tm_scale_bar(position = c('right', 'bottom'), width = 0.3)

Gi_race <- localG(CD_race_sf$X..of.population.non.white, spatial_weights, zero.policy = TRUE)
CD_Gi_race <- cbind(nyc_community_districts, as.matrix(Gi_race))  
names(CD_Gi_race)[5] <- 'Getis Ord Gi*'

tm_shape(CD_Gi_race) +
  tm_fill("Getis Ord Gi*",
          palette = '-RdBu',
          style = 'pretty') +
  tm_borders(alpha = 0.2) +
  tm_layout(main.title = "Hot and cold spots of racial distribution",
            main.title.position = 'centre',
            frame = FALSE) +
  tm_compass(position = c('right', 'bottom')) +
  tm_scale_bar(position = c('right', 'bottom'), width = 0.3)
  
Gi_income <- localG(CD_income_sf$Median_Income, spatial_weights, zero.policy = TRUE)
CD_Gi_income <- cbind(nyc_community_districts, as.matrix(Gi_income))
names(CD_Gi_income)[5] <- 'Getis Ord Gi*'

tm_shape(CD_Gi_income) +
  tm_fill("Getis Ord Gi*",
          palette = '-RdBu',
          style = 'pretty') +
  tm_borders(alpha = 0.2) +
  tm_layout(main.title = "Hot and cold spots of income distribution",
            main.title.position = 'centre',
            frame = FALSE) +
  tm_compass(position = c('right', 'bottom')) +
  tm_scale_bar(position = c('right', 'bottom'), width = 0.3)

