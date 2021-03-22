#TAKE A LOOK AT SAMPLE DATA FROM Y.Li
#WRITTEN BY SR, SPRING 2021

# Load everything ---------------------------------------------------------
library(tidyverse)
theme_set(theme_bw()) #Good for maps
library(sf)
library(ggmap)
library(ggpubr)

#Water data
wDat <- read.csv('./data/sample_waterData.csv') %>% mutate(Date=as.Date(Date,'%Y-%m-%d')) %>% 
  st_as_sf(coords=c('lon','lat')) %>% st_set_crs(4326)

#Get locations from wDat to join onto spatial data
locIndex <- wDat %>% select(YEID,geometry) %>% distinct()

#Spectral data
sDat <- read.csv('./data/sample_spectralData.csv') %>% 
  left_join(locIndex,by='YEID') %>% #Join in spatial info
  mutate(date_img=as.Date(date_img,'%Y-%m-%d'))

#Base map
spDomain <- wDat %>% st_union() %>% st_convex_hull() %>% st_buffer(dist=1) %>% st_bbox() #Spatial domain
names(spDomain) <- c('left','bottom','right','top')
basemap <- get_map(location=spDomain) #Get from Google Earth

# Take a look at data -----------------------------------------------------

#Water data

#only surface water
surfWDat <- wDat %>% group_by(YEID) %>% mutate(depthOrd=order(Depth)) %>% ungroup() %>% 
  filter(depthOrd==1) %>% select(-depthOrd)

# ggmap(basemap)+ #Split up by date
ggplot()+
  geom_sf(data=surfWDat, aes(col=DO),inherit.aes = FALSE)+
  facet_wrap(~Date)+
  scale_colour_distiller(type='seq',palette = "YlOrBr")

p1 <- ggplot(surfWDat)+geom_sf(aes(col=DO))+
  scale_colour_distiller(type='div',palette = "YlOrBr")
p2 <- ggplot(surfWDat)+geom_sf(aes(col=Salin))+
  scale_colour_distiller(type='div',palette = "YlOrBr")
p3 <- ggplot(surfWDat)+geom_sf(aes(col=Temp))+
  scale_colour_distiller(type='div',palette = "YlOrBr")
ggarrange(p1,p2,p3,ncol=1)

#Basic bivariate correlation plots
panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
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

surfWDat %>% st_drop_geometry() %>% select(DO,Temp,Salin) %>% pairs(upper.panel=panel.cor,diag.panel=panel.hist)

#Time series plots - may have to restrict analysis to specific time "chunks"
surfWDat %>% st_drop_geometry() %>% pivot_longer(c(DO,Temp,Salin)) %>% 
  ggplot(aes(x=Date,y=value))+geom_point()+
  facet_wrap(~name,ncol=1)

#Next: look at spectral data and decide how these could be put together

