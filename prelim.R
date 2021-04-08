#TAKE A LOOK AT SAMPLE DATA FROM Y.Li
#WRITTEN BY SR, SPRING 2021

# Load everything ---------------------------------------------------------
library(tidyverse)
theme_set(theme_bw()) #Good for maps
library(sf)
library(ggmap)
library(ggpubr)
library(animation)

setwd("~/Documents/hypoxiaMapping")

#Water data
wDat <- read.csv('./data/sample_waterData.csv') %>% mutate(Date=as.Date(Date,'%Y-%m-%d')) %>% 
  # mutate(YEID=ifelse(YEID=='2014_149','2014_003',YEID)) %>% #Set 2014_149 as 2014_003 (see below)
  group_by(YEID) %>% mutate(maxDepth=max(Depth,Depth_dem)) %>% mutate(propDepth=Depth/maxDepth) %>%
  select(-Depth_dem) %>% relocate(YEID:Date,Depth,maxDepth,propDepth,DO,Temp:lon) %>% 
  st_as_sf(coords=c('lon','lat')) %>% st_set_crs(4326)

# # 2014_003, 2014_149 are different times at the same location
# wDat %>% mutate(x=st_coordinates(.)[,1],y=st_coordinates(.)[,2]) %>% unite(loc,x,y) %>% st_drop_geometry() %>% select(YEID,loc) %>%
#   distinct() %>% group_by(loc) %>% mutate(n=n()) %>% filter(n>1)

#Get locations from wDat to join onto spatial data
locIndex <- wDat %>% select(YEID,geometry) %>% distinct()

#Spectral data
sDat <- read.csv('./data/sample_spectralData.csv') %>% 
  left_join(locIndex,by='YEID') %>% #Join in spatial info
  st_as_sf() %>% 
  mutate(date_img=as.Date(date_img,'%Y-%m-%d'))

#Not sure what the other channels are, so removing them for now, and stripping out data-less days
sDat2 <- sDat %>% select(YEID:poc,geometry) %>% 
  mutate(noData=is.na(chlor_a)&is.na(nflh)&is.na(poc)) %>% 
  group_by(date_img) %>% mutate(n=n(),nNoDat=sum(noData)) %>% #All have 205 locations
  filter(nNoDat<(n*0.5)) %>% #Remove days where >50% of data are missing
  select(-noData:-nNoDat) %>% ungroup() %>% 
  mutate(nflh=ifelse(nflh<0,0,ifelse(nflh>1,1,nflh))) %>% #Limits between 0 and 1
  mutate(poc=ifelse(nflh<0,1,poc)) #Changes negative poc to 1 (really low)

#Base map
spDomain <- wDat %>% st_union() %>% st_convex_hull() %>% st_buffer(dist=1) %>% st_bbox() #Spatial domain
names(spDomain) <- c('left','bottom','right','top')
basemap <- get_map(location=spDomain) #Get from Google Earth

#choose only surface water, using order of depth measurements
surfWDat <- wDat %>% group_by(YEID) %>% mutate(depthOrd=order(Depth)) %>% ungroup() %>% 
  filter(depthOrd==1) %>% select(-depthOrd)

#choose only bottom water
bottomWDat <- wDat %>% group_by(YEID) %>% mutate(depthOrd=order(Depth,decreasing=TRUE)) %>% ungroup() %>% 
  filter(depthOrd==1) %>% select(-depthOrd)

# Helper functions -----------

midcut <- function(dat,b,...){ #Function to cut data and assign label as midpoint (not boundaries)
  x <- cut(dat,breaks=b,...) # cut the data into bins...
  if(length(b)==1){ #If single break value
    s <- seq(min(dat,na.rm=TRUE),max(dat,na.rm=TRUE),length.out=b+1) #Sequence of boundary points
  } else { #If multiple break values
    s <- b; b <- length(s)
  }
  levels(x) <- round(rowMeans(cbind(s,lag(s)))[2:c(b+1)],1) #Average with lag, then assign as new levels
  return(x)
}

dateCut <- function(d,b,include.lowest=TRUE){ #Cut date object d into b discrete chunks, and return as date object
  #include.lowest = should date on left boundary be included in category 1 (otherwise returns NA)
  if(class(b)=='character' & length(b)>1){ 
    b <- as.numeric(as.Date(b)) 
  } else if(!(length(b)==1 & class(b)=='numeric')) {
    stop('b must be single number, or character vector of date breaks in YYYY-MM-DD format')
  }
  cdates <- midcut(as.numeric(d),b,include.lowest=include.lowest) #Convert date to numeric, and chop into pieces
  cdates <- round(as.numeric(as.character(cdates))) #Turn into integers (using midpoint between breaks)
  cdates <- as.Date(cdates,origin="1970-01-01") #Turn into dates
  return(cdates)
}

#Turn sf geometry (x/y coordinates) into columns
#   removeGeom: should geometry be dropped?
#   epsgCode: should CRS be transformed (useful for lat/lon -> UTM)
geom2cols <- function(d,x=lon,y=lat,removeGeom=TRUE,epsg=NA){
  if(!is.na(epsg)) d <- st_transform(d,epsg) #Transform to new CRS
  d <- d %>% mutate({{x}}:=st_coordinates(.)[,1],{{y}}:=st_coordinates(.)[,2]) #Make new columns from coordinates
  if(removeGeom) d <- st_drop_geometry(d) #Drop geometry
  return(d)
}

# Take a look at water data ---------------------

#Overall
p1 <- ggmap(basemap)+geom_sf(data=surfWDat,aes(col=DO),inherit.aes = FALSE)+
  scale_colour_distiller(type='div',palette = "YlOrBr")
p2 <- ggmap(basemap)+geom_sf(data=surfWDat,aes(col=Salin),inherit.aes = FALSE)+ #Salin
  scale_colour_distiller(type='div',palette = "YlOrBr")
p3 <- ggmap(basemap)+geom_sf(data=surfWDat,aes(col=Temp),inherit.aes = FALSE)+ #Temp
  scale_colour_distiller(type='div',palette = "YlOrBr")
p <- ggarrange(p1,p2,p3,ncol=1)
ggsave('./figures/prelimFigs/wDat_overall_surf.png',p,width=8,height=8)

p1 <- ggmap(basemap)+geom_sf(data=bottomWDat,aes(col=DO),inherit.aes = FALSE)+
  scale_colour_distiller(type='div',palette = "YlOrBr")
p2 <- ggmap(basemap)+geom_sf(data=bottomWDat,aes(col=Salin),inherit.aes = FALSE)+ #Salin
  scale_colour_distiller(type='div',palette = "YlOrBr")
p3 <- ggmap(basemap)+geom_sf(data=bottomWDat,aes(col=Temp),inherit.aes = FALSE)+ #Temp
  scale_colour_distiller(type='div',palette = "YlOrBr")
p <- ggarrange(p1,p2,p3,ncol=1)
ggsave('./figures/prelimFigs/wDat_overall_bottom.png',p,width=8,height=8)

#Split up maps by date
p <- ggplot()+
  geom_sf(data=surfWDat, aes(col=DO),inherit.aes = FALSE)+
  facet_wrap(~Date)+
  scale_colour_distiller(type='seq',palette = "YlOrBr")
ggsave('./figures/prelimFigs/wDat_DO_date.png',p,width=8,height=6)

#Time series plots - may have to restrict analysis to specific time "chunks"
p <- surfWDat %>% st_drop_geometry() %>% pivot_longer(c(DO:Salin)) %>% 
  ggplot(aes(x=Date,y=value))+geom_point()+
  facet_wrap(~name,ncol=1)
ggsave('./figures/prelimFigs/wDat_overall_ts.png',p,width=8,height=8)

#Basic bivariate correlation plots
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "grey", ...)
}
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * abs(r),col=ifelse(r>0,'black','red'))
}

png(file = './figures/prelimFigs/wDat_overall_cor.png',width=6,height=6,units='in',res=150)
surfWDat %>% st_drop_geometry() %>% select(DO:Salin) %>% pairs(upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()

#How do water values change with depth?
# Salin increases with depth (mixing at surface near rivers), temp decreases (vertical mixing)
# DO decreases with depth, but highly variable in general
p <- wDat %>% st_drop_geometry() %>% 
  pivot_longer(c(DO:Salin)) %>% mutate(fDate=dateCut(Date,6)) %>% 
  ggplot()+  geom_point(aes(x=value,y=Depth),alpha=0.3)+ #Uses absolute depth
  facet_grid(fDate~name,scales='free_x')+
  scale_y_reverse()+labs(y='Depth',x='Measurement')
ggsave('./figures/prelimFigs/wDat_depth_overall.png',p,width=8,height=10)

wDat %>% st_drop_geometry() %>% 
  pivot_longer(c(DO:Salin)) %>% mutate(fDate=dateCut(Date,6)) %>% 
  ggplot()+  geom_point(aes(x=value,y=Depth/Depth_dem),alpha=0.3)+ #Uses proportional depth
  facet_grid(fDate~name,scales='free_x')+
  scale_y_reverse()+labs(y='Proportion Depth',x='Measurement')

wDat %>% group_by(YEID) %>% mutate(mDepth=max(Depth,Depth_dem)) %>% data.frame()

wDat %>% group_by(YEID) %>% slice(1)


# Take a look at spectral data ----------------------------
sDat
head(sDat)

p <- sDat %>% select(YEID:poc) %>% pivot_longer(c(chlor_a:poc)) %>% filter(!is.na(value)) %>% 
  ggplot()+geom_point(aes(x=date_img,y=value))+facet_wrap(~name,ncol=1,scales='free_y')+
  geom_vline(xintercept=range(wDat$Date),col='red',linetype='dashed') #Range of water data
ggsave('./figures/prelimFigs/sDat_overall_ts.png',p,width=8,height=8)


#Make animation of spectral data across time range
gen_anim <- function() { #Generate set of frames
  sDat3 <- sDat2 %>% mutate(across(chlor_a:poc,~as.vector(scale(.x)))) #%>% pivot_longer(chlor_a:poc)
  mapFun <- function(tau) { #Function to make single frame
    s <- sDat3 %>%  filter(date_img == tau) #subset data to time step tau
      # %>% filter(!is.na(value))
    f <- function(dat,channel,lims){ #Function to make standard figure
      ggmap(basemap) + geom_sf(data=s,aes(geometry=geometry,colour = {{channel}},size={{channel}}),inherit.aes = FALSE) +
        guides(size=FALSE)+
        scale_size(limits=lims)+
        # scale_colour_gradient(limits=lims)
        scale_colour_distiller(limits=lims,type='div',palette = "YlOrBr",na.value=NA)
    }
    p1 <- f(s,chlor_a,range(sDat3$chlor_a,na.rm=TRUE))
    p2 <- f(s,nflh,range(sDat3$nflh,na.rm=TRUE))
    p3 <- f(s,poc,range(sDat3$poc,na.rm=TRUE))
    p <- ggarrange(p1,p2,p3,ncol=1)
    p <- annotate_figure(p,top = text_grob(as.character(tau),size=10))
    return(p)
  }
  mapFun(tau=sDat2$date_img[1])

  for(t in 1:length(unique(sDat2$date_img))){  # for each time point
    print(mapFun(tau=unique(sDat2$date_img)[t]))   # plot data at this time point
  }
}

setwd("~/Documents/hypoxiaMapping/figures/prelimFigs")
saveGIF(gen_anim(),movie.name='spectral_anim.gif',interval = 0.5,ani.width=600,ani.height=900)
setwd("~/Documents/hypoxiaMapping")

#Hard to see distinct patterns using animations, but it looks like there are "pulses" here and there

#Hovmoller plots for each channel
p1 <- sDat2 %>% mutate(lon=st_coordinates(.)[,1]) %>% 
  mutate(lon=midcut(lon,20)) %>% st_drop_geometry() %>% 
  mutate(across(chlor_a:poc,~scale(.x))) %>% 
  group_by(date_img,lon) %>% summarize(across(c(chlor_a:poc),mean,na.rm=TRUE)) %>% ungroup() %>% 
  mutate(lon=as.numeric(as.character(lon))) %>% 
  pivot_longer(chlor_a:poc) %>% 
  ggplot(aes(x=lon,y=date_img))+geom_tile(aes(fill=value))+facet_wrap(~name)+
  scale_fill_distiller(type='div',palette = "YlOrBr",na.value=NA)+
  theme(panel.background = element_rect(fill = 'grey'))+
  labs(x='Longitude',y='Date')+theme(legend.position = 'bottom')

p2 <- sDat2 %>% mutate(lat=st_coordinates(.)[,2]) %>% 
  mutate(lat=midcut(lat,20)) %>% st_drop_geometry() %>% 
  mutate(across(chlor_a:poc,~scale(.x))) %>% 
  group_by(date_img,lat) %>% summarize(across(c(chlor_a:poc),mean,na.rm=TRUE)) %>% ungroup() %>% 
  mutate(lat=as.numeric(as.character(lat))) %>% 
  pivot_longer(chlor_a:poc) %>% 
  ggplot(aes(x=lat,y=date_img))+geom_tile(aes(fill=value))+facet_wrap(~name)+
  scale_fill_distiller(type='div',palette = "YlOrBr",na.value=NA)+
  theme(panel.background = element_rect(fill = 'grey'))+
  labs(x='Latitude',y='Date')+theme(legend.position = 'bottom')
p <- ggarrange(p1,p2,ncol=2,common.legend=TRUE,legend='right')
ggsave('./figures/prelimFigs/sDat_hovmoller.png',p,width=12,height=8)

# sDat2 %>% mutate(lon=st_coordinates(.)[,1]) %>% mutate(lon=cut(lon,breaks=20)) %>% st_drop_geometry() %>% 
#   group_by(date_img,lon) %>% summarize(across(c(chlor_a:poc),mean,na.rm=TRUE)) %>% ungroup() %>% pull(lon) %>% 
#   levels(.) %>% str()

#Pretty low correlation between classes. Do these actually mean the same thing? Ask YL or LN how these are typically used.
png(file = './figures/prelimFigs/sDat_overall_cor.png',width=6,height=6,units='in',res=150)
sDat2 %>% select(chlor_a:poc) %>% na.omit() %>% st_drop_geometry() %>% pairs(upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()

# Fit model of chlor_a ----------------------------------------------------
library(mgcv)

chlorDat <- sDat2 %>% geom2cols(E,N,epsg=3401) %>% 
  mutate(sE=as.vector(scale(E)),sN=as.vector(scale(N))) %>% 
  mutate(doy=format(date_img,format='%j'),doy=as.vector(scale(as.numeric(doy)))) %>% 
  mutate(logChlor=log(chlor_a)) %>% 
  filter(!is.na(chlor_a))

# library(parallel) #No significant difference in fitting times using bam
# detectCores()
# cl <- makeCluster(12)
# stopCluster(cl) 

chlorMod <- gam(logChlor~te(sN,sE,doy,k=9),data=chlorDat)  

summary(chlorMod)
par(mfrow=c(2,2)); gam.check(chlorMod); par(mfrow=c(1,1))
plot(chlorMod,scheme=3,too.far=0.1)

chlorDat$resid <- resid(chlorMod)

ggplot(chlorDat)+geom_point(aes(x=E,y=N,col=resid,size=abs(resid)))+
  scale_colour_distiller(type='div',palette = "YlOrBr")


