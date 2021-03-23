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
  st_as_sf(coords=c('lon','lat')) %>% st_set_crs(4326)

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

# Take a look at data -----------------------------------------------------

#Water data

#only surface water
surfWDat <- wDat %>% group_by(YEID) %>% mutate(depthOrd=order(Depth)) %>% ungroup() %>% 
  filter(depthOrd==1) %>% select(-depthOrd)

#Overall
p1 <- ggmap(basemap)+geom_sf(data=surfWDat,aes(col=DO),inherit.aes = FALSE)+
  scale_colour_distiller(type='div',palette = "YlOrBr")
p2 <- ggmap(basemap)+geom_sf(data=surfWDat,aes(col=Salin),inherit.aes = FALSE)+ #Salin
  scale_colour_distiller(type='div',palette = "YlOrBr")
p3 <- ggmap(basemap)+geom_sf(data=surfWDat,aes(col=Temp),inherit.aes = FALSE)+ #Temp
  scale_colour_distiller(type='div',palette = "YlOrBr")
p <- ggarrange(p1,p2,p3,ncol=1)
ggsave('./figures/prelimFigs/wDat_overall.png',p,width=8,height=8)

#Split up maps by date
p <- ggplot()+
  geom_sf(data=surfWDat, aes(col=DO),inherit.aes = FALSE)+
  facet_wrap(~Date)+
  scale_colour_distiller(type='seq',palette = "YlOrBr")
ggsave('./figures/prelimFigs/wDat_DO_date.png',p,width=8,height=6)

#Time series plots - may have to restrict analysis to specific time "chunks"
p <- surfWDat %>% st_drop_geometry() %>% pivot_longer(c(DO,Temp,Salin)) %>% 
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
surfWDat %>% st_drop_geometry() %>% select(DO,Temp,Salin) %>% pairs(upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()

# Spectral data
sDat
head(sDat)

p <- sDat %>% select(YEID:poc) %>% pivot_longer(c(chlor_a:poc)) %>% filter(!is.na(value)) %>% 
  ggplot()+geom_point(aes(x=date_img,y=value))+facet_wrap(~name,ncol=1,scales='free_y')+
  geom_vline(xintercept=range(wDat$Date),col='red',linetype='dashed') #Range of water data
ggsave('./figures/prelimFigs/sDat_overall_ts.png',p,width=8,height=8)



#Make animation of spectral data across time range 
gen_anim <- function() {
  mapFun <- function(tau) { #Function to make single frame
    s <- sDat2 %>%  filter(date_img == tau)   # subset data to time step tau
    f <- function(dat,channel,lims){ #Function to make standard figure
      ggmap(basemap) + geom_sf(data=s,aes(geometry=geometry,colour = {{channel}},size={{channel}}),inherit.aes = FALSE) + 
        guides(size=FALSE)+
        scale_size(limits=lims)+
        scale_colour_gradient(limits=lims)
    }
    p1 <- f(s,chlor_a,range(sDat2$chlor_a,na.rm=TRUE))
    p2 <- f(s,nflh,range(sDat2$nflh,na.rm=TRUE))
    p3 <- f(s,poc,range(sDat2$poc,na.rm=TRUE))
    p <- ggarrange(p1,p2,p3,ncol=1)
    p <- annotate_figure(p,top = text_grob(as.character(tau),size=10))
    return(p)
  }
  # test <- mapFun(tau)
  
  for(t in 1:length(unique(sDat2$date_img))){  # for each time point
    print(mapFun(tau=unique(sDat2$date_img)[t]))   # plot data at this time point
  }
}

setwd("~/Documents/hypoxiaMapping/figures/prelimFigs")
saveGIF(gen_anim(),movie.name='spectral_anim.gif',interval = 0.5,ani.width=600,ani.height=900)
setwd("~/Documents/hypoxiaMapping")

#Hard to see distinct patterns using animations, but it looks like there are "pulses" here and there

#Hovmoller plots for each channel

midcut<-function(dat,b){ #Function to cut data and assign label as midpoint (not boundaries)
  x <- cut(dat,breaks=b) # cut the data into bins...
  s <- seq(min(dat,na.rm=TRUE),max(dat,na.rm=TRUE),length.out=b+1) #Sequence of boundary points
  levels(x) <- round(rowMeans(cbind(s,lag(s)))[2:c(b+1)],1) #Average with lag, then assign as new levels
  return(x)
}

sDat2 %>% mutate(lon=st_coordinates(.)[,1]) %>% mutate(lon=midcut(lon,20)) %>% st_drop_geometry() %>% 
  group_by(date_img,lon) %>% summarize(across(c(chlor_a:poc),mean,na.rm=TRUE)) %>% ungroup() %>% 
  mutate(lon=as.numeric(as.character(lon))) %>% 
  pivot_longer(chlor_a:poc) %>% 
  ggplot(aes(x=lon,y=date_img))+geom_tile(aes(fill=value))+facet_wrap(~name)

sDat2 %>% mutate(lon=st_coordinates(.)[,1]) %>% mutate(lon=cut(lon,breaks=20)) %>% st_drop_geometry() %>% 
  group_by(date_img,lon) %>% summarize(across(c(chlor_a:poc),mean,na.rm=TRUE)) %>% ungroup() %>% pull(lon) %>% 
  levels(.) %>% str()

sDat2 %>% mutate(lon=st_coordinates(.)[,1]) %>% pull(lon) %>% min()



#Pretty low correlation between classes. Do these actually mean the same thing? Ask YL or LN how these are typically used.
png(file = './figures/prelimFigs/sDat_overall_cor.png',width=6,height=6,units='in',res=150)
sDat2 %>% select(chlor_a:poc) %>% na.omit() %>% pairs(upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()




