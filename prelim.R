#TAKE A LOOK AT SAMPLE DATA FROM Y.Li
#WRITTEN BY SR, SPRING 2021

# Load everything ---------------------------------------------------------
library(tidyverse)
theme_set(theme_bw()) #Good for maps
library(sf)
library(ggmap)
library(ggpubr)
library(animation)
library(mgcv)

setwd("~/Documents/hypoxiaMapping")

source('helperFunctions.R')

#Water data
wDat <- read.csv('./data/sample_waterData.csv') %>% mutate(Date=as.Date(Date,'%Y-%m-%d')) %>% mutate(doy=as.numeric(format(Date,format='%j'))) %>% 
  # mutate(YEID=ifelse(YEID=='2014_149','2014_003',YEID)) %>% #Set 2014_149 as 2014_003 (see below)
  group_by(YEID) %>% mutate(maxDepth=max(Depth,Depth_dem)) %>% mutate(propDepth=Depth/maxDepth) %>%
  select(-Depth_dem) %>% relocate(YEID:Date,doy,Depth,maxDepth,propDepth,DO,Temp:lon) %>% ungroup() %>% 
  st_as_sf(coords=c('lon','lat')) %>% st_set_crs(4326) %>% 
  geom2cols(E,N,removeGeom=FALSE,epsg=3401) %>% #Louisiana offshore 
  mutate(sE=(E-mean(E))/1000,sN=(N-mean(N))/1000) %>% #Center E/N and convert to km
  st_transform(4326)

# #All YEIDs are from a single day
# wDat %>% st_drop_geometry() %>% group_by(YEID) %>% summarize(nDays=length(unique(Date))) %>% filter(nDays!=1)

# # 2014_003, 2014_149 are different times at the same location
# wDat %>% mutate(x=st_coordinates(.)[,1],y=st_coordinates(.)[,2]) %>% unite(loc,x,y) %>% st_drop_geometry() %>% select(YEID,loc) %>%
#   distinct() %>% group_by(loc) %>% mutate(n=n()) %>% filter(n>1)

#Get locations from wDat to join onto spatial data
locIndex <- wDat %>% select(YEID,Date,doy,E:sN,geometry) %>% distinct() 
  
#Spectral data
sDat <- read.csv('./data/sample_spectralData.csv') %>% 
  left_join(locIndex,by='YEID') %>% #Join in spatial info
  select(-Date,-doy) %>% 
  st_as_sf() %>% 
  mutate(date_img=as.Date(date_img,'%Y-%m-%d'),doy=as.numeric(format(date_img,'%j'))) %>% 
  mutate(nflh=ifelse(nflh<0,0,ifelse(nflh>1,1,nflh))) %>% #Limits nflh to between 0 and 1
  mutate(poc=ifelse(nflh<0,1,poc)) #Changes negative poc to 1 (really low)

#Not sure what the other channels are, so removing them for now, and stripping out data-less days
sDat2 <- sDat %>% select(YEID,date_img,doy,chlor_a:poc,E:geometry) %>% 
  mutate(noData=is.na(chlor_a)&is.na(nflh)&is.na(poc)) %>% #Strip days where all data are missing
  group_by(date_img) %>% mutate(n=n(),nNoDat=sum(noData)) %>% #All have 205 locations
  filter(nNoDat<(n*0.5)) %>% #Remove days where >50% of data are missing
  select(-noData:-nNoDat) %>% ungroup() 

#Base map
spDomain <- wDat %>% st_union() %>% st_convex_hull() %>% st_buffer(dist=0.25) %>% st_bbox() #Spatial domain
names(spDomain) <- c('left','bottom','right','top')
basemap <- get_map(location=spDomain) #Get from Google Earth

#choose only surface water, using order of depth measurements
surfWDat <- wDat %>% group_by(YEID) %>% mutate(depthOrd=order(Depth)) %>% ungroup() %>% 
  filter(depthOrd==1) %>% select(-depthOrd)

#choose only bottom water
bottomWDat <- wDat %>% group_by(YEID) %>% mutate(depthOrd=order(Depth,decreasing=TRUE)) %>% ungroup() %>% 
  filter(depthOrd==1) %>% select(-depthOrd)

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

#Same thing, but with proportion depth
wDat %>% st_drop_geometry() %>% 
  pivot_longer(c(DO:Salin)) %>% mutate(fDate=dateCut(Date,6)) %>% 
  ggplot()+  geom_point(aes(x=value,y=propDepth),alpha=0.3)+ #Uses proportional depth
  facet_grid(fDate~name,scales='free_x')+
  scale_y_reverse()+labs(y='Proportion Depth',x='Measurement')

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
      ggmap(basemap) + geom_sf(data=dat,aes(geometry=geometry,colour = {{channel}},size={{channel}}),inherit.aes = FALSE) +
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
(p <- ggarrange(p1,p2,ncol=2,common.legend=TRUE,legend='right'))
ggsave('./figures/prelimFigs/sDat_hovmoller.png',p,width=12,height=8)

# sDat2 %>% mutate(lon=st_coordinates(.)[,1]) %>% mutate(lon=cut(lon,breaks=20)) %>% st_drop_geometry() %>% 
#   group_by(date_img,lon) %>% summarize(across(c(chlor_a:poc),mean,na.rm=TRUE)) %>% ungroup() %>% pull(lon) %>% 
#   levels(.) %>% str()

#Pretty low correlation between classes. Do these actually mean the same thing? Ask YL or LN how these are typically used.
png(file = './figures/prelimFigs/sDat_overall_cor.png',width=6,height=6,units='in',res=150)
sDat %>% select(chlor_a:poc) %>% na.omit() %>% st_drop_geometry() %>% pairs(upper.panel=panel.cor,diag.panel=panel.hist)
dev.off()

#Almost no chlor_a measurements on days when DO measurments were taken. Need some way of filling in the gaps
bottomWDat %>% left_join(mutate(st_drop_geometry(sDat),doy=doy),by=c('YEID','doy')) %>% 
  filter(!is.na(chlor_a)) %>% 
  ggplot(aes(chlor_a,DO))+geom_point()

# Fit GAM of chlor_a ----------------------------------------------------

#Using entire dataset. Fit is somewhat better if using reduced dataset (sDat2), but going with this for now.
#Would fit be better using reduced dataset because of 
chlorDat <- sDat %>% select(YEID:chlor_a,E:geometry) %>% 
  mutate(doy=format(date_img,format='%j'),doy=as.numeric(doy)) %>% 
  mutate(logChlor=log(chlor_a)) %>% 
  filter(!is.na(chlor_a)) #Strips out missing data

# library(parallel) #No significant difference in fitting times using bam
# detectCores()
# cl <- makeCluster(12)
# stopCluster(cl) 

#Uses 50 bases for N,E and 10 for time. Tensor product = 500 terms
chlorMod <- gam(logChlor~te(sN,sE,doy,bs=c('tp','tp'),k=c(50,10),d=c(2,1)),data=chlorDat)
# chlorMod <- gam(chlor_a~te(sN,sE,doy,bs=c('tp','tp'),k=c(50,10),d=c(2,1)),data=chlorDat,family=Gamma) #Takes a ton of time
# chlorMod2 <- gam(chlor_a~te(sN,sE,doy,bs=c('tp','tp'),k=c(30,5),d=c(2,1)),data=chlorDat,family=Gamma())
# chlorMod3 <- gam(chlor_a~te(sN,sE,doy,bs=c('tp','tp'),k=c(30,5),d=c(2,1)),data=chlorDat,family=inverse.gaussian()) 
# par(mfrow=c(2,2)); gam.check(chlorMod2); abline(0,1,col='red'); par(mfrow=c(1,1))
# par(mfrow=c(2,2)); gam.check(chlorMod3); abline(0,1,col='red'); par(mfrow=c(1,1))

summary(chlorMod)
par(mfrow=c(2,2)); gam.check(chlorMod); abline(0,1,col='red'); par(mfrow=c(1,1))
plot(chlorMod,scheme=3,too.far=0.1)
chlorDat$resid <- resid(chlorMod) #Residuals

#Prediction grid
predChlor <- with(chlorDat,expand.grid(doy=seq(min(doy),max(doy),length.out=9),loc=unique(paste(sE,sN,sep='_')))) %>% 
  separate(loc,c('sE','sN'),sep='_',convert=TRUE) %>% 
  mutate(pred=predict(chlorMod,newdata=.),se=predict(chlorMod,newdata=.,se.fit=TRUE)$se.fit) %>% 
  mutate(date=as.Date(paste('2020',round(doy),sep='-'),format='%Y-%j'))

p1 <- predChlor %>% ggplot()+geom_point(aes(x=sE,y=sN,col=pred))+
  facet_wrap(~date,nrow=1)+
  scale_colour_distiller(type='div',palette = "Greens",direction=1)+
  labs(x='E',y='N',col='Predicted\nlog(Chlor)')+
  theme(legend.position='bottom')

p2 <- predChlor %>% ggplot()+geom_point(aes(x=sE,y=sN,col=se))+
  facet_wrap(~date,nrow=1)+
  scale_colour_distiller(type='div',palette = "YlOrRd",direction=1) +
  labs(x='E',y='N',col='SE')+
  theme(legend.position='bottom')

(p <- ggarrange(p1,p2,nrow=2))
ggsave('./figures/prelimFigs/chlorMod.png',p,width=11,height=6)

#Plot of residuals over space for 12 time slices
chlorDat %>% mutate(cDate=dateCut(date_img,12)) %>%
  ggplot()+geom_sf(aes(geometry=geometry,col=resid,size=abs(resid)),alpha=0.5)+
  facet_wrap(~cDate)+ scale_colour_distiller(type='div',palette = "RdBu")

#Variogram of residuals

library(gstat) #Spatial-only variogram

#Overall spatial autocorrelation
select(chlorDat,YEID,date_img,resid) %>% 
  st_transform(3401) %>% as('Spatial') %>% 
  variogram(resid~1, data=.,width=2000) %>% 
  ggplot()+geom_line(aes(x=dist,y=gamma))+
  labs(x='Distance',y='Semivariance')

chlorDat %>% filter(date_img==first(date_img)) %>% 
  select(YEID,date_img,resid) %>% 
  st_transform(3401) %>% as('Spatial') %>% 
  variogram(resid~1, data=.,width=2000) %>% 
  ggplot()+geom_line(aes(x=dist,y=gamma))+
  labs(x='Distance',y='Semivariance')

chlorDat %>% filter(date_img==last(date_img)) %>% 
  select(YEID,date_img,resid) %>% 
  st_transform(3401) %>% as('Spatial') %>% 
  variogram(resid~1, data=.,width=2000) %>% 
  ggplot()+geom_line(aes(x=dist,y=gamma))+
  labs(x='Distance',y='Semivariance')

#Spatial autocorr (Moran's I) at each day - strong spatial autocorrelation
m <- sapply(sort(unique(chlorDat$date_img)),function(i){
  temp <- chlorDat %>% filter(date_img==i)
  
  if(nrow(temp)<=5) return(NA)
  
  d <- temp %>% st_distance()
  units(d) <- NULL
  d <- 1000*d/max(d) #Change to km
  d <- 1/d #inverse
  diag(d) <- 0
  r <- temp %>% pull(resid) #Residuals
  if(sum(!is.finite(d))>0){ #IF NA/Inf (sites at same location)
    notFinite <- which(apply(d,2,function(x) as.logical(sum(!is.finite(x)))))[-1]
    d <- d[-notFinite,-notFinite] #Remove rows with NA/Inf
    r <- r[-notFinite]
  }
  if(!is.finite(Moran.I(r,d)$p.value)) stop('NaN value')
  
  return(Moran.I(r,d)$p.value)
})
m <- m[!is.na(m)] #Remove NAs
sum(m<(0.05/length(m)))/length(m) #73% of days show evidence of spatial autocorrelation

#Temporal autocorrelation (Durbin-Watson) at each location - uses order rather than day
n <- sapply(unique(chlorDat$YEID),function(i){
  temp <- chlorDat %>% filter(YEID==i) %>% dwtest(resid~1,data=.)
  return(temp$p.value)
})
sum(n<(0.05/length(n)))/length(n) #Not much

# #ACF of residuals. Creates ts() object with NAs in missing days
# temp <- chlorDat %>% filter(YEID==last(YEID)) %>% st_drop_geometry() %>% mutate(doy=1+doy-min(doy))
# r <- rep(NA,max(temp$doy))
# r[temp$doy] <- temp$resid
# a <- acf(r,na.action = na.pass,type='partial',lag.max=30)
# 
# r <- ts(r,start=1,end=length(r))
# ar(r,na.action='na.pass')
# Box.test(r)

#Empirical ST variogram of residuals
a <- Sys.time()
vv <- chlorDat %>% slice_sample(n=1000) %>% 
  st_drop_geometry() %>% 
  stConstruct(x=.,space=c('E','N'),time='date_img') %>% 
  variogramST(resid~1, #Takes about 2 mins if using 1000 data points, width=2000, tlags= 100, cores =10
            data=.,
            width=10000,
            tlags=seq(0,30,by=1),cores = 10,tunit='days')
Sys.time()-a

vv %>% 
  ggplot(aes(x=spacelag/1000,y=timelag,fill=gamma))+geom_raster() +
  labs(x='Distance (km)',y='Time (days)',fill='Semi-\nvariance')



# #Try larger tensor product (10x10x10).  What has 500 more coefficients gotten us?
# chlorMod2 <- gam(logChlor~te(sN,sE,doy,k=10),data=chlorDat)  #Better AIC, but residual checks are about the same.
# summary(chlorMod)
# summary(chlorMod2)
# anova(chlorMod,chlorMod2,test='Chisq') #Slightly worse CV (0.347 vs 0.350)
# AIC(chlorMod,chlorMod2) #AIC is worse
# chlorDat$resid2 <- resid(chlorMod2)
# 
# #Plot of residuals over space for 5 time slices. Looks like both models are doing poorly at the same places
# pivot_longer(chlorDat,resid:resid2) %>%
#   mutate(cDate=dateCut(date_img,6)) %>%
#   ggplot()+geom_point(aes(x=E,y=N,col=value,size=abs(value)),alpha=0.5)+
#   facet_grid(cDate~name)+
#   scale_colour_distiller(type='div',palette = "RdBu")
# 
# #Mod1 is worse around the edges of the distribution, while mod2 is worse at the centre
# chlorDat %>% mutate(resDiff=abs(resid)-abs(resid2)) %>%  #Difference of absolute residuals
#   mutate(cDate=dateCut(date_img,12)) %>%
#   ggplot()+geom_point(aes(x=E,y=N,col=resDiff,size=abs(resDiff)),alpha=0.5)+
#   facet_wrap(~cDate)+
#   scale_colour_distiller(type='div',palette = "RdBu") #Red = model1 worse, Blue = model2 worse
# 
# #Looks about the same in time
# pivot_longer(chlorDat,resid:resid2) %>%
#   ggplot(aes(x=date_img,y=value))+geom_point()+
#   facet_wrap(~name,ncol=1)

# Predict DO using chlor_a ------------------------------------------------

#Simple model using lagged chlor_a to predict DO for both 

#Get predictions at all locations through the entire season
predChlor <- with(chlorDat,expand.grid(YEID=unique(YEID),doy=min(doy):max(doy))) %>% 
  left_join(select(locIndex,-doy,-Date),by='YEID') %>% #Join in spatial info
  mutate(pred=predict(chlorMod,newdata=.),se=predict(chlorMod,newdata=.,se.fit=TRUE)$se.fit) %>% 
  mutate(date=as.Date(paste('2020',round(doy),sep='-'),format='%Y-%j'))

lags <- -40:40 #Try -40 to 40 day lags (negative = days after, positive = days before)
modList <- lapply(lags,function(i){
  #Copies of water data
  sWatTemp <- surfWDat
  bWatTemp <- bottomWDat 
  chooseThis <- predChlor$YEID==sWatTemp$YEID & predChlor$doy==(sWatTemp$doy-i) #Location of lagged chlor predictions
  sWatTemp$lagChlor <- predChlor$pred[chooseThis]
  sWatTemp$lagDay <- predChlor$doy[chooseThis]
  bWatTemp$lagChlor <- predChlor$pred[chooseThis]
  bWatTemp$lagDay <- predChlor$doy[chooseThis]
  sMod <- lm(DO~lagChlor,data=sWatTemp)
  bMod <- lm(DO~lagChlor,data=bWatTemp)
  return(list(surface=sMod,bottom=bMod))
})

p1 <- data.frame(lag=lags,surface=sapply(modList,function(i) sum(resid(i$surface)^2)),
           bottom=sapply(modList,function(i) sum(resid(i$bottom)^2))) %>% 
  pivot_longer(surface:bottom) %>% 
  ggplot()+geom_line(aes(x=lag,y=value))+
  geom_vline(xintercept = 0,linetype='dashed')+facet_wrap(~name)+
  labs(x='Time lag',y='Mean Squared Error')

p2 <- data.frame(lag=lags,surface=sapply(modList,function(i) summary(i$surface)$r.squared),
           bottom=sapply(modList,function(i) summary(i$bottom)$r.squared)) %>% 
  pivot_longer(surface:bottom) %>% 
  ggplot()+geom_line(aes(x=lag,y=value))+
  geom_vline(xintercept = 0,linetype='dashed')+facet_wrap(~name)+
  labs(x='Time lag',y='R-squared')

ggarrange(p1,p2,ncol=1)

#Take a look at model on day -11
par(mfrow=c(2,1))
dayLag <- which.min(sapply(modList,function(i) sum(resid(i$bottom)^2)))
bestMod <- modList[[dayLag]]$bottom #Bottom water model
with(bestMod$model,plot(lagChlor,DO,xlab=paste0('Chlor_a at day ',lags[dayLag]),ylab='Bottom DO'))
abline(bestMod) #Fairly good relationship (not linear, but OK for now)

dayLag <- which.min(sapply(modList,function(i) sum(resid(i$surface)^2)))
bestMod <- modList[[dayLag]]$surface #Surface water model
with(bestMod$model,plot(lagChlor,DO,xlab=paste0('Chlor_a at day ',lags[dayLag]),ylab='Surface DO'))
abline(bestMod) #No real relationship
par(mfrow=c(1,1))

# NOTES:
# Looks like surface water DO is more poorly predicted by chlor_a. However, lower MSE, so probably just less variable overall
# Strangely, future chlor_a is a slightly better predictor than past chlor_a. What is going on here?
# Hypotheses:
# 1) Just a weird dataset, and won't show up in future years
# 2) Some thing is actually causing a DO drop, followed by an increase in chlor_a


#Trying functional regression using lagged predictions
NdayLag <- 30 #30 days in past
NdayForward <- 0 #20 days in future
dayLags <- -NdayForward:NdayLag

#Matrix to store chlor_a predictions for past 0:30 days
predChlorMat <- matrix(NA,nrow=nrow(bottomWDat),ncol=length(dayLags),
                       dimnames=list(bottomWDat$YEID,gsub('-','m',paste0('lag',dayLags))))

for(i in 1:nrow(bottomWDat)){
  getDays <- bottomWDat$doy[i]:bottomWDat$doy[i]-dayLags
  predChlorMat[i,] <- predChlor$pred[predChlor$YEID == bottomWDat$YEID[i] & predChlor$doy %in% getDays]
}

#Data for functional regression
fdat <- list(DO_bottom=bottomWDat$DO,DO_surf=surfWDat$DO,
          dayMat=outer(rep(1,nrow(bottomWDat)),dayLags),
          chlorMat=predChlorMat,
          doy=bottomWDat$doy,sE=bottomWDat$sE,sN=bottomWDat$sN,
          maxDepth=bottomWDat$maxDepth)

#Fit FDA models - looks slightly better, but still only gets R2 of 0.3, so not much better than simple lag model
bWatMod <- gam(DO_bottom ~ s(dayMat,by=chlorMat), data=fdat) #Bottom water
summary(bWatMod)
par(mfrow=c(2,2)); gam.check(bWatMod); abline(0,1,col='red'); par(mfrow=c(1,1)) #Not too bad
plot(bWatMod,scheme=1,pages=1); abline(h=0,col='red')

data.frame(pred=predict(bWatMod),actual=fdat$DO_bottom) %>% 
  ggplot()+geom_point(aes(x=pred,y=actual))+
  geom_abline(intercept = 0, slope = 1)+
  labs(x='Predicted Bottom DO',y='Actual Bottom DO')

p1 <- bottomWDat %>% mutate(resid=resid(bWatMod)) %>% st_jitter(factor=0.01) %>% 
  ggplot()+geom_sf(aes(geometry=geometry,col=resid,size=abs(resid)),alpha=0.5)+
  scale_colour_distiller(type='div',palette = "RdBu")+theme(legend.position='bottom')+
  guides(size=FALSE)
p2 <- bottomWDat %>% mutate(resid=resid(bWatMod)) %>% ggplot()+geom_point(aes(x=doy,y=resid))+geom_hline(yintercept = 0)
ggarrange(p1,p2,ncol=1)

sWatMod <- gam(DO_surf ~ s(dayMat,by=chlorMat), data=fdat) #Surface water DO
summary(sWatMod)
par(mfrow=c(2,2)); gam.check(sWatMod); abline(0,1,col='red'); par(mfrow=c(1,1))
plot(sWatMod,pages=1,scheme=1); abline(h=0,col='red') 

p1 <- surfWDat %>% mutate(resid=resid(sWatMod)) %>% st_jitter(factor=0.01) %>% 
  ggplot()+geom_sf(aes(geometry=geometry,col=resid,size=abs(resid)),alpha=0.5)+
  scale_colour_distiller(type='div',palette = "RdBu")+theme(legend.position='bottom')+
  guides(size=FALSE)
p2 <- surfWDat %>% mutate(resid=resid(sWatMod)) %>% ggplot()+geom_point(aes(x=doy,y=resid))+geom_hline(yintercept = 0)
ggarrange(p1,p2,ncol=1)

#To try:
#PCA of chlor_a, nflh, poc + other channels

