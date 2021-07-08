#CREATE PREDICTED BOTTOM DO MAPS BASED ON SURFACE MEASUREMENTS
#WRITTEN BY SR, SUMMER 2021

#Load everything -------------------------------------------------

library(tidyverse)
theme_set(theme_bw()) #Good for maps
library(sf)
library(raster)
source('helperFunctions.R')

storage <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/ATdata"
newfolder <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/ATdata/combined"

#Get file paths
files <- data.frame(paths=dir(storage,recursive=TRUE,full.names = TRUE)) %>% 
  mutate(file=dir(storage,recursive=TRUE)) %>%
  filter(!grepl('combined',paths)) %>% 
  mutate(file=gsub('.*\\/','',file)) %>% mutate(doy=as.numeric(gsub('\\.tif','',gsub('[AT]2014','',file)))) %>%
  dplyr::select(-file) %>% 
  mutate(platform=ifelse(grepl('terr',paths),'terra','aqua')) %>% 
  mutate(doy=ifelse(platform=='terra',doy-1,doy)) %>% #Terra starts 1 day earlier
  pivot_wider(names_from='platform',values_from = paths) %>%
  mutate(combined=gsub('\\/media.*\\/[AT]','',aqua)) %>% 
  mutate(combined=gsub('2014','2014_',combined)) %>% 
  mutate(combined=paste0(newfolder,'/',combined)) %>% 
  mutate(exists=file.exists(combined))

# a <- Sys.time()
# for(i in 1:nrow(files)){
#   if(!files$exists[i]){
#     aquaDat <- brick(files$aqua[i])
#     terraDat <- brick(files$terra[i])
#     
#     newDat <- lapply(1:dim(aquaDat)[3],function(j) calc(stack(aquaDat[[j]],terraDat[[j]]),fun=mean,na.rm=TRUE))
#     newDat <- brick(newDat)
#     names(newDat) <- names(aquaDat)
#     writeRaster(newDat,filename=files$combined[i],format='GTiff')
#     rm(newDat)
#     gc()    
#   }
#   print(paste0('Finished file ',i,' of ',nrow(files)))
# }
# b <- Sys.time()
# b-a #Takes about an hour. Also causes problems because of background file storage within raster library

brick(files$aqua[1])

sDat <- brick(files$combined[1])
