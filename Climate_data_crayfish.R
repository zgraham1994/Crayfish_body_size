
#Crayfish Climate Date Extraction with Raster
##################
#### packages ####
##################
library(dplyr)
library(purrr)
library(readr)  
library(magrittr)
library(rgbif)
library(taxize)
library(rgdal)
library(raster)
library(sp)
library(phytools)
library(geiger)
library(caper)
library(ggplot2)
library(ggtree)
library(ggnewscale)

data<-read.csv("ContainedCentroid.csv", row.names=1, stringsAsFactors = TRUE)

####################################
##### CLIMATE DATA ACQUISITION #####
####################################

#we just want species name and location, reduce dataset to that
cols.to.keep <- c("species","decimalLatitude","decimalLongitude")
Clade1_occ <- data[, cols.to.keep]
#plot data to get an idea of location data
plot(Clade1_occ$decimalLongitude, Clade1_occ$decimalLatitude, pch=".")
#remove occurrences that have 0 for lat or long
m1 <- subset(Clade1_occ, Clade1_occ$decimalLatitude != 0.00)
m1 <- subset(m1, m1$decimalLongitude != 0.00)

Clade1_occ$decimalLongitude <- as.numeric(as.character(Clade1_occ$decimalLongitude))
Clade1_occ$decimalLatitude <- as.numeric(as.character(Clade1_occ$decimalLatitude))
#write
#write.csv(Clade1_occ, file = "", row.names = F)

# designating the coordinates for raster::extract
coordinates(Clade1_occ) <- ~ decimalLongitude + decimalLatitude
#setwd("C:/Users/Desktop/wc2.1_2.5m_tavg")
#setwd("C:/Users/Zack/Desktop/wc2.1_2.5m_tavg")
#I had to pul all of these files form the WC zip file onto my desktop, set the WD to desktop
#and then it compiled them together
tmean_files <- list.files(path = ".", pattern = "wc2.1_")

# Import using the purrr::map function to apply the raster::stack function
tmean_rasters <- raster::stack(tmean_files)

# adding climate data to species list
# the numbers coorespond to the lat/lon columns of the data frame
BIO1 <- extract(tmean_rasters, Clade1_occ@coords)
tmean.avg <- rowMeans(BIO1)

Clade1_df <- cbind.data.frame(Clade1_occ, tmean.avg)
Clade1_df <- na.omit(Clade1_df)
write.csv(Clade1_df, "updated.csv", row.names = F)
