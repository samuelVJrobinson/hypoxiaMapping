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
  require(sf); require(dplyr)
  if(is_grouped_df(d)){ #If data are grouped
    warning('Data are grouped. Will ungroup data before converting geometry column.')
    d <- ungroup(d)
  }
  if(!is.na(epsg)) d <- st_transform(d,epsg) #Transform to new CRS
  d <- d %>% mutate({{x}}:=st_coordinates(.)[,1],{{y}}:=st_coordinates(.)[,2]) #Make new columns from coordinates
  if(removeGeom) d <- st_drop_geometry(d) #Drop geometry
  return(d)
}
