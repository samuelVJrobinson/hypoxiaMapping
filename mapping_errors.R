#Look into why the mapping predictions are so far off

# Load data -----------------

library(tidyverse)
theme_set(theme_bw()) #Good for maps
library(sf)
library(raster)
library(mgcv)
library(parallel)
source('helperFunctions.R')

load('./data/all2014_2.Rdata')

sDat_model <- sDat #Model data

rm(bottomWDat,locIndex,sDat,surfWDat,wDat)

storage <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/ATdata"
newfolder <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/ATdata/combined"
shpfileFolder <- "/media/rsamuel/Storage/geoData/Rasters/hypoxiaMapping2021/shapefiles"

coastBuff <- st_read(paste0(shpfileFolder,"/region_coast250kmBuffer.shp")) #250 km buffer zone away from coast
# st_read(paste0(shpfileFolder,"/region_rectangle.shp")) %>%
# st_read(paste0(shpfileFolder,"/TWAP_rivers.shp")) %>%
# ggplot()+geom_sf()

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

sDat <- lapply(files$combined,function(datPath){
  cNames <- names(brick(files$aqua[1])) #Channel names
  sD <- brick(datPath) #Read data
  names(sD) <- cNames #Add names
  sD2 <- rasterToPoints(sD) #Convert to point dataframe
  nMissing <- apply(sD2[,!grepl('(x|y)',colnames(sD2))],1,function(x) sum(is.na(x))) #Proportion of missing values in each row (cell)
  sD2 <- data.frame(doy=strsplit(x = datPath,split = c("(\\_|\\.)"))[[1]][2],sD2[nMissing<5,]) #Keep points with 4 or more, and add date
  return(sD2)}) #Takes about 10 seconds
names(sDat) <- sapply(sDat,function(x) x$doy[1])

# #Preview day 152
# sDat[[1]] %>% #Chlor_a on day 152
#   # filter(x>(-86.75),x<(-86.25),y<30.2,y>29.8) %>%
#   ggplot(aes(x=x,y=y,fill=chlor_a))+
#   geom_raster()

#Decided to omit non-complete cells - lots of completely missing cells already, so this doesn't hurt that much
sDat2 <- do.call('rbind',sDat) %>% na.omit() %>% #Combine into single DF and remove NAs
  mutate(doy=as.numeric(doy)) #Convert to numeric

#I think this is causing the problem:  
sDat2 <- sDat2 %>% 
  # mutate(nflh=rescale(nflh,1e-5,(1-1e-5))) %>% #Rescales nflh to between 0 and 1
  mutate(nflh=ifelse(nflh>(1-1e-5),(1-1e-5),nflh)) %>% #Set upper limit of nflh just below 1
  mutate(across(chlor_a:sst,~ifelse(.x<0,lwrLimits[names(lwrLimits)==cur_column()]*0.95,.x))) #Rescales negative values be above 0.95*minimum positive value
summary(sDat2)

#Only observations
obs <- sDat2 %>% dplyr::select(-doy:-y) %>% as.matrix() %>% log()
# (obs[c(1),]-pca1$center)/pca1$scale %*% pca1$rotation #Works with a single row

#Calculate PCs 1-6
# ((sDatMat_imputed$completeObs[1,]-pca1$center)/pca1$scale) %*% pca1$rotation 
# ((obs[1,]-pca1$center)/pca1$scale) %*% pca1$rotation
PCs <- ((obs-outer(rep(1,nrow(obs)),pca1$center))/outer(rep(1,nrow(obs)),pca1$scale)) %*% pca1$rotation
PCs <- PCs[,1:6]
sDat2 <- cbind(sDat2,PCs) #Combine PCs with sDat2
# sDat2 %>% st_as_sf(coords=c("x","y")) %>% filter(doy==153) %>% ggplot()+geom_sf(aes(col=PC1))

#Add coordinate system
sDat2 <- sDat2 %>% st_as_sf(coords=c('x','y')) %>% st_set_crs(4326) %>%
  geom2cols(E,N,removeGeom=FALSE,epsg=3401) %>% #Louisiana offshore
  mutate(sE=(E-mean(unique(E)))/1000,sN=(N-mean(unique(N)))/1000) %>% #Center E/N and convert to km
  mutate(across(E:N,~round(.x))) %>% 
  st_transform(4326) %>% unite(loc,E,N,remove = FALSE) %>% 
  mutate(loc=as.numeric(factor(loc)))

withinBuff <- sDat2 %>% st_intersects(.,coastBuff) %>% sapply(.,function(x) length(x)>0) #Points in sDat2 that are outside of the buffer
sDat2 <- sDat2 %>% filter(withinBuff) #Filter out points outside of buffer
locLookup <- sDat2 %>% dplyr::select(loc:geometry) %>% unique() #Lookup table for locations
# ggplot(locLookup)+geom_sf()

#Keeping geometry for now
sDat2 <- sDat2 %>% mutate(doy=as.Date(paste0('2014-',doy),format='%Y-%j')) %>% 
  rename(date_img=doy)# %>% st_drop_geometry()
# sDat_model <- sDat_model %>% st_drop_geometry()

# Figures -------------
head(sDat2)
head(sDat_model)

#Raw values - scaling of nflh might be causing problems
temp1 <- sDat2 %>% st_drop_geometry() %>% dplyr::select(date_img:sst) %>% 
  pivot_longer(cols = -date_img) %>% 
  filter(!is.na(value)) %>% mutate(type='mapping')
temp2 <- sDat_model %>% st_drop_geometry() %>% dplyr::select(date_img:sst) %>% 
  pivot_longer(cols = -date_img) %>% 
  filter(!is.na(value)) %>% mutate(type='modeling')

(p1 <- bind_rows(temp1,temp2) %>% 
  group_by(date_img,name,type) %>% 
  summarise(med=median(value),max=quantile(value,0.9),min=quantile(value,0.1)) %>% 
  ggplot(aes(x=date_img,y=med,col=type))+
  geom_pointrange(alpha=0.5,size=0.5,aes(ymax=max,ymin=min))+
  facet_wrap(~name,scales='free_y')+
  # scale_y_log10()+
  labs(x='Date',y='Value (min-median-max)',title='Raw spectral data',col=NULL)+
  theme(legend.position = 'bottom'))

temp1 <- sDat2 %>% st_drop_geometry() %>% dplyr::select(date_img,contains('PC')) %>% 
  pivot_longer(cols = -date_img) %>% 
  filter(!is.na(value)) %>% mutate(type='mapping')
temp2 <- sDat_model %>% st_drop_geometry() %>% dplyr::select(date_img,contains('PC')) %>% 
  pivot_longer(cols = -date_img) %>% 
  filter(!is.na(value)) %>% mutate(type='modeling')

(p2 <- bind_rows(temp1,temp2) %>% 
  group_by(date_img,name,type) %>% 
  summarise(med=median(value),max=quantile(value,0.9),min=quantile(value,0.1)) %>% 
  ggplot(aes(x=date_img,y=med,col=type))+
  geom_pointrange(alpha=0.5,size=0.5,aes(ymax=max,ymin=min))+
  facet_wrap(~name,scales='free_y')+
  # scale_y_log10()+
  labs(x='Date',y='Value (min-median-max)',title='Principle components',col=NULL)+
  theme(legend.position = 'bottom'))

ggsave(p1,filename = './figures/errors_raw2.png',width=16,height=10)
ggsave(p2,filename = './figures/errors_PC2.png',width=16,height=10)

#Regression plots

head(sDat2) #Mapping data
head(sDat_model) #Modeling data

#Restrict model data to date range of mapping data
temp_modDat <- sDat_model %>% filter(date_img>=min(sDat2$date_img) & date_img<=max(sDat2$date_img)) 

modDatBuff <- temp_modDat %>% st_union() %>% st_convex_hull() %>% st_buffer(dist = 0.5) #0.5 degree buffer around modeling data

temp_locLookup <- slice(locLookup,unlist(st_intersects(modDatBuff,locLookup))) #locLookup locations inside buffer

yeid_lookup <- temp_modDat %>% dplyr::select(YEID) %>% distinct() #Distinct locations for modeling data

yeid_lookup$loc <- temp_locLookup$loc[apply(st_distance(yeid_lookup,temp_locLookup),1,which.min)] #Nearest location in locLookup

temp_modDat <- left_join(temp_modDat,st_drop_geometry(yeid_lookup),by='YEID') %>% mutate(ID=paste(date_img,loc,sep='_'))  #Join in location index, create ID column

temp_mapDat <- sDat2 %>% st_drop_geometry() %>% mutate(ID=paste(date_img,loc,sep='_')) %>% dplyr::select(ID,chlor_a:sst)

bothDat <- temp_modDat %>% dplyr::select(chlor_a:sst,ID) %>% 
  left_join(temp_mapDat,by='ID',suffix=c('.modDat','.mapDat')) %>% 
  st_drop_geometry() %>% group_by(ID) %>% summarise(across(everything(),mean)) %>% ungroup() %>% 
  pivot_longer(cols=-ID) %>% 
  separate(col=name,into=c('variable','dataset'),sep='\\.') %>% arrange(ID,dataset,variable) %>% 
  pivot_wider(names_from=dataset,id_cols=ID:variable,values_from=value) %>% 
  na.omit()

(p1 <- bothDat %>% ggplot(aes(x=mapDat,y=modDat))+
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype='dashed',col='red')+
  facet_wrap(~variable,scales='free')+
  labs(x='Mapping data (new)',y='Modeling data (original)'))

ggsave(p1,filename = './figures/compare_data2.png',width=16,height=10)
