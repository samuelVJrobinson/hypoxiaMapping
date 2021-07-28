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

# #Average between aqua and terra
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
# b-a #Takes about 10 mins

#Combined channels are missing their names, but this isn't that big of a deal. 
cNames <- names(brick(files$aqua[1])) #Channel names
datPath <- files$combined[1] #First combined stack
# datPath <- files$aqua[1] #Aqua stack

sDat <- brick(datPath)
names(sDat) <- cNames
sDat2 <- rasterToPoints(sDat) #Convert to points
nMissing <- apply(sDat2[,!grepl('(x|y)',colnames(sDat2))],1,function(x) sum(is.na(x))) #Proportion of missing values in each row (cell)
sDat2 <- data.frame(sDat2[nMissing<5,]) #Keep points with 4 or more 

sDat2 %>% 
  filter(x>(-86.75),x<(-86.25),y<30.2,y>29.8) %>% #Weird sampling going on; way too many points per cell
  ggplot(aes(x=x,y=y,fill=chlor_a))+
  geom_raster()

#No simple way to do this using raster. Resample using tidyverse syntax?

# ext <- extent(c(-86.65,-86.35,29.8,30.1))
# sDatSmall <- crop(sDat[[1]],ext)
# 
# plot(sDatSmall)
# extract(sDatSmall,ext)
# 
# # convert to SPDF
# gridDat <- as(sDatSmall,'SpatialPolygonsDataFrame')
# spplot(gridDat)
# 
# sDatSmall2 <- aggregate(sDat,fact=c(10,9),fun=max)
# 
# plot(sDatSmall)

l <- rep(c(1,2,5,3),c(3,4,2,6))
l2 <- rep(c(3,NA,3.5,4,NA,5),c(5,2,3,5,5,5))

blockSize <- function(x){ #Gets size of blocks in a vector of NAs/reals
  xNA <- is.na(x)

  x1 <- x[!xNA]-lag(x[!xNA],1)
  x1[1] <- 1
  x1 <- x1!=0
  x2 <- as.numeric(!xNA)
  x2[!xNA] <- cumsum(x1)
  
  x3 <- x2 - lag(x2,1)
  x3[1] <- 1
  x3 <- x3!=0
  unname(table(cumsum(x3)))
}

debugonce(blockSize)
blockSize(l)
blockSize(l2)


apply(as.matrix(sDat[[1]]),1,blockSize)
# apply(as.matrix(sDat[[1]]),2,blockSize)

blockGroup <- function(x){
  xNA <- is.na(x)
  x1 <- x[!xNA]-lag(x[!xNA],1)
  x1[1] <- 1
  x1[abs(x1)<.Machine$double.eps] <- 0
  x1 <- x1!=0
  cumsum(x1)
  x[!xNA] <- cumsum(x1)
  return(x)
}
# 
# blockGroup(l)
# blockGroup(l2)

sDat3 <- sDat2 %>% 
  arrange(x,y) %>% mutate(across(.cols=chlor_a:sst,~blockGroup(.x),.names = "{.col}_xgroup")) %>%
  unite(col='x_group',contains("_xgroup"),sep='_') %>% mutate(x_group=as.numeric(factor(x_group))) %>%
  arrange(y,x) %>% mutate(across(.cols=chlor_a:sst,~blockGroup(.x),.names = "{.col}_ygroup")) %>%
  unite(col='y_group',contains("_ygroup"),sep='_') %>% mutate(y_group=as.numeric(factor(y_group))) %>%
  unite(col='cell_group',contains("_group"),sep='_') %>% mutate(cell_group=as.numeric(factor(cell_group))) %>% 
  group_by(cell_group) %>% summarize(across(everything(),~mean(.x)))

sDat3 %>% filter(x>(-86.75),x<(-86.25),y<30.2,y>29.8) %>% #Weird sampling going on; way too many points per cell
  ggplot(aes(x=x,y=y,col=chlor_a))+geom_point(size=0.5) 

sDat2 %>% 
  # filter(x>(-86.75),x<(-86.25),y<30.2,y>29.8) %>%
  filter(x>(-86.4),x<(-86.3),y<30,y>29.9) %>%
  arrange(y,x) %>% mutate(across(.cols=chlor_a:sst,~blockGroup(.x),.names = "{.col}_xgroup")) %>%
  unite(col='x_group',contains("_xgroup"),sep='_') %>% mutate(x_group=as.numeric(factor(x_group))) %>%
  arrange(x,y) %>% mutate(across(.cols=chlor_a:sst,~blockGroup(.x),.names = "{.col}_ygroup")) %>%
  unite(col='y_group',contains("_ygroup"),sep='_') %>% mutate(y_group=as.numeric(factor(y_group))) %>%
  ggplot(aes(y=y))+
  # geom_text(aes(x=x,label=y_group),alpha=0.5,col='red')+
  geom_text(aes(x=x+0.001,label=x_group),col='black')
  # geom_point(aes(col=y_group),size=1.5) 
