### MARTIN ET AL. 2025
## Life history traits influence environmental impacts on spatial population synchrony in European birds and butterflies
## Ecological Monologues
## Code supplement: Calculating synchrony in the environment (temperature and precipitation)
## Date modified: 29-01-2025

country = "replace with countryname"
taxa = "replace with taxa"
method = "AvgGrid"

## Load packages
library(dplyr)
library(reshape2)
library(tidyr)
library(geosphere)
library(sf)
library(mapsf)
library(rgdal)
library(mvtnorm)
library(RColorBrewer)
library(data.table)
library(scales)


## Define  functions for estimating synchrony ####
GaussSyncE <- function(rho.0, rho.inf, x, lambda) (rho.0-rho.inf)*exp(-(x^2)/(2*lambda^2)) + rho.inf #Sjekk parantesen i siste term
ExpSyncE <- function(rho.0, rho.inf, x, lambda) (rho.0-rho.inf)*exp(-x/lambda) + rho.inf # Sjekk parantesen i siste term

GaussSyncF <- function(par, x, distmat, startVal=FALSE){
  # x : standardisert
  
  rho.0 <- exp(par[1])/(1+exp(par[1]))
  rho.inf <- rho.0*exp(par[2])/(1+exp(par[2]))          
  lambda <- exp(par[3])
  if (startVal){
    cat(round(c(rho.0, rho.inf, lambda)),3)
  }
  logLL <- 0
  for (i in 1:nrow(x)){
    xx <- x[i,]
    sel <- !is.na(xx)
    distmatsel <- distmat[sel,sel]
    xx <- xx[sel]
    mu <- rep(0, length(xx))
    sigma <- (rho.0-rho.inf)*exp(-(distmatsel^2)/(2*lambda^2)) + rho.inf
    diag(sigma) <- 1
    logLL <- logLL + -sum(dmvnorm(xx, mu, sigma,log = TRUE))
  }
  logLL
}
GaussSyncFboot <- function(par, x, distmat, startVal=FALSE){
  # x : standardisert
  
  rho.0 <- exp(par[1])/(1+exp(par[1]))
  rho.inf <- rho.0*exp(par[2])/(1+exp(par[2]))          
  lambda <- exp(par[3])
  if (startVal){
    cat(round(c(rho.0, rho.inf, lambda)),3)
  }
  logLL <- 0
  for (i in 1:nrow(x)){
    xx <- x[i,]
    sel <- !is.na(xx)
    distmatsel <- distmat[sel,sel]
    xx <- xx[sel]
    mu <- rep(0, length(xx))
    sigma <- (rho.0-rho.inf)*exp(-(distmatsel^2)/(2*lambda^2)) + rho.inf
    #diag(sigma) <- 1
    logLL <- logLL + -sum(dmvnorm(xx, mu, sigma,log = TRUE))
  }
  logLL
}

GaussSyncFFix <- function(par, x, distmat, startVal=FALSE){
  # x : standardisert
  
  logLL <- 0
  for (i in 1:nrow(x)){
    xx <- x[i,]
    sel <- !is.na(xx)
    distmatsel <- distmat[sel,sel]
    xx <- xx[sel]
    mu <- rep(0, length(xx))
    sigma <- matrix(par, nr=sum(sel), nc=sum(sel))
    diag(sigma) <- 1
    logLL <- logLL + -sum(dmvnorm(xx, mu, sigma,log = TRUE))
  }
  logLL
}
ExpSyncF <- function(par, x, distmat, startVal=FALSE){
  # x : standardisert
  
  rho.0 <- exp(par[1])/(1+exp(par[1]))
  rho.inf <- rho.0*exp(par[2])/(1+exp(par[2]))          
  lambda <- exp(par[3])
  if (startVal){
    cat(round(c(rho.0, rho.inf, lambda)),3)
  }
  logLL <- 0
  for (i in 1:nrow(x)){
    xx <- x[i,]
    sel <- !is.na(xx)
    distmatsel <- distmat[sel,sel]
    xx <- xx[sel]
    mu <- rep(0, length(xx))
    sigma <- (rho.0-rho.inf)*exp(-distmatsel/lambda) + rho.inf
    diag(sigma) <- 1
    logLL <- logLL + -sum(dmvnorm(xx, mu, sigma,log = TRUE))
  }
  logLL
}


## Begin analysis ####

counts <- read.csv("count file location") ## Raw counts per survey area
sites <-  read.csv("survey location information") 
names(sites)[1]<-"SITE.CODE"

## Merge site locations to count indices file
CountsSites <- merge(counts, sites, by='SITE.CODE')
UniqueSurveyLocs <- CountsSites %>% distinct(SITE.CODE, .keep_all=TRUE)

## Define projections for sf objects. Layers are ETRS and BNG: want to transform to UTM33
ETRSproj = "+init=epsg:4258"
BNGproj = "+init=epsg:27700" ## Region boundaries in British national grid projection
UTM33proj <- "+proj=utm +zone=33 +datum=WGS84"


## Create grid boundaries
Country <- sf::st_read("\\NUTS_LVL_3\\Europe_NUTSLVL3.shp") ## Import all European country boundaries to use to create grid
## NUTS levels are accessed from here: https://ec.europa.eu/eurostat/web/gisco/geodata/statistical-units/territorial-units-statistics
Country_NUTS3 <- Country[Country$LEVL_CODE == "3",] ## Select the THIRD NUTS level
COUNTRY_FLAT <- st_transform(Country_NUTS3, crs=3035)
tempgrid <- st_make_grid(COUNTRY_FLAT, cellsize = c(100000,100000), square = FALSE, crs = 3035) %>%   ##Create an overlay 100 by 100km grid hexagons
  st_sf() %>% 
  mutate(id_100km = 1:nrow(.))## Create hexagonal grid over europe study countries


## Make the site count data sf objects, transform to correct crs (lonlatproj)
CountData <- st_as_sf(x=CountsSites, coords=c("Easting", "Northing"), crs=BNGproj)
DataUTM <- st_transform(CountData, crs=3035) ## Transform combined sf object into ETRS projection.
SiteLocs <- st_as_sf(x=UniqueSurveyLocs, coords=c("Easting", "Northing"), crs=BNGproj)
SiteLocsUTM <- st_transform(SiteLocs, crs=3035) ## Transform combined sf object into ETRS projection.


## Spatial overlay to assign points to regions/counties
overlaypoints <- st_intersects(DataUTM, tempgrid)
overlaypoints_2 <- unlist(lapply(overlaypoints, function(x) ifelse(length(x)==0, NA, as.vector(tempgrid$id_100km)[x])))
DataUTM[,"region"] <- overlaypoints_2


# Aggregate the counts by the regions/counties (NAME)
DataAgg <- aggregate(SITE.INDEX~region+SPECIES+YEAR, data=DataUTM, 
                     FUN=function(x) sum(x, na.rm=TRUE))
DataAggWide <- spread(DataAgg, key=YEAR, value=SITE.INDEX) ## Transform the dataframe to a wide format


## Calculate averages
table <- tapply(DataUTM$SITE.CODE, list(DataUTM$YEAR, DataUTM$region), function(x) length(unique(x)))
widetable <- as.data.frame(t(table))
widetable <- cbind(rownames(widetable),widetable)
rownames(widetable) <- NULL
colnames(widetable) <- c(names(widetable)) #to not write all the column names
colnames(widetable)[1] <- "region" 
names(widetable)
YEAR <- "YEAR"
numbersites <- "numbersites"
## Gather all the columns (years) which have abundance data
gathercols <- c("1973","1974","1975","1976","1977","1978","1979","1980","1981","1982","1983","1984","1985","1986","1987",
                "1988","1989","1990", "1991","1992", "1993", "1994","1995","1996","1997","1998","1999","2000","2001", 
                "2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013","2014","2015", "2016","2017", "2018", "2019", "2020")
longtable<- tidyr::gather(widetable, YEAR, numbersites, gathercols)
DataAgg$region <- as.character(DataAgg$region)
longtable$YEAR <- as.character(longtable$YEAR)
DataAgg$YEAR <- as.character(DataAgg$YEAR)

DataAggAVG <- DataAgg %>%
  left_join(longtable) 

## Calculate average below
DataAggAVG$average <- with(DataAggAVG, SITE.INDEX/numbersites)

DataAggAVGWide <- subset(DataAggAVG, select = -c(SITE.INDEX, numbersites))
DataAggAVGWide<- reshape(DataAggAVGWide, idvar = c("region","SPECIES"), timevar = "YEAR", direction = "wide")
names(DataAggAVGWide)[3:ncol(DataAggAVGWide)] <- c("1973","1974","1975","1976","1977","1978","1979","1980","1981","1982","1983","1984","1985","1986","1987",
                                                   "1988","1989","1990", "1991","1992", "1993", "1994","1995","1996","1997","1998","1999","2000","2001", 
                                                   "2002","2003","2004","2005","2006",
                                                   "2007","2008","2009","2010","2011","2012","2013","2014","2015", "2016","2017", "2018", "2019", "2020")
names(DataAggAVGWide)[names(DataAggAVGWide) == "1990"] <- "X1990"
names(DataAggAVGWide)[names(DataAggAVGWide) == "1991"] <- "X1991"
names(DataAggAVGWide)[names(DataAggAVGWide) == "1992"] <- "X1992"
names(DataAggAVGWide)[names(DataAggAVGWide) == "1993"] <- "X1993"
names(DataAggAVGWide)[names(DataAggAVGWide) == "1994"] <- "X1994"
names(DataAggAVGWide)[names(DataAggAVGWide) == "1995"] <- "X1995"
names(DataAggAVGWide)[names(DataAggAVGWide) == "1996"] <- "X1996"
names(DataAggAVGWide)[names(DataAggAVGWide) == "1997"] <- "X1997"
names(DataAggAVGWide)[names(DataAggAVGWide) == "1998"] <- "X1998"
names(DataAggAVGWide)[names(DataAggAVGWide) == "1999"] <- "X1999"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2000"] <- "X2000"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2000"] <- "X2001"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2002"] <- "X2002"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2003"] <- "X2003"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2004"] <- "X2004"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2005"] <- "X2005"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2006"] <- "X2006"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2007"] <- "X2007"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2008"] <- "X2008"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2009"] <- "X2009"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2010"] <- "X2010"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2011"] <- "X2011"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2012"] <- "X2012"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2013"] <- "X2013"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2014"] <- "X2014"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2015"] <- "X2015"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2016"] <- "X2016"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2017"] <- "X2017"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2018"] <- "X2018"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2019"] <- "X2019"
names(DataAggAVGWide)[names(DataAggAVGWide) == "2020"] <- "X2020"

## Cleaning Data: Grid
## Unique regions post-aggregation to be used as new populations:
aggregated_lans <- unique(DataUTM$region)

## Remove locations where, per species, there are fewer than 50 individuals observed across all years.
DataAggAVGWide$total <- rowSums(DataAggAVGWide[3:ncol(DataAggAVGWide)], na.rm=TRUE)
DataAggWide <- DataAggAVGWide[DataAggAVGWide[,ncol(DataAggAVGWide)] >= 50, ]
DataAggWideSubset <- subset(DataAggWide, select = -c(total))

## Removing all species per population where there are more than 5 years with missing data (NAs).
DataAggWideSubset <- DataAggWideSubset[rowSums(is.na(DataAggWideSubset)) <= 5, ]


## Set initial species list for population dynamics plotting
Species <- unique(DataAggWideSubset$SPECIES)
unique.species.subset <- Species
unique.species.subset_numbered<- 1:length(unique.species.subset)
names(unique.species.subset_numbered) <- unique.species.subset
names(unique.species.subset_numbered)


## Setting up synchrony calculations
## Use only regions with more than 1 sampling point
npoints <- table(DataAggWideSubset$region)
subset <- npoints[npoints>1]
LocsSubset <- DataUTM[DataUTM$region%in%names(subset),]
DataAggSubset <- DataAggWideSubset[DataAggWideSubset$region%in%names(subset),]

## Create coordinate and distance matrix
Coords <- aggregate(st_coordinates(LocsSubset), 
                    by=list(region=data.frame(LocsSubset)[,"region"]), FUN=mean)
distmat <- dist(Coords[,c("X", "Y")], upper=TRUE)/1000


SyncData <- DataAggSubset[DataAggSubset$SPECIES%in%Species,] # check for the exclusion of poor sampling years
SyncDataList <- split(SyncData, f=list(
  SyncData$SPECIES, SyncData$region), drop=TRUE)
SyncDataList <- lapply(SyncDataList, function(x) 
  reshape(x, direction="long", varying=list(3:ncol(SyncData)), sep=""))
SyncDataList <- lapply(SyncDataList, function(x) 
  cbind.data.frame(x, log.r1=append(NA, diff(log(x$X))))) 
SyncDataList <- lapply(SyncDataList, function(x) 
  cbind.data.frame(x, log.r=ifelse(is.na(x$log.r1)|is.nan(x$log.r1)|is.infinite(x$log.r1), NA, x$log.r1)))
SyncDataList <- do.call(rbind.data.frame, SyncDataList)
SyncDataList <- split(SyncDataList, f=SyncDataList$SPECIES)

# make wide format
SyncDatalogr <- lapply(SyncDataList, function(x) 
  reshape(x[,c("time", "region", "log.r")], direction="wide", 
          idvar="time", timevar="region"))
SyncDatalogr <- lapply(SyncDatalogr, function(x) 
  scale(x[apply(x, 1, function(y) sum(!is.na(y))>1),-c(1)]))
SyncLists <- list(log.r=SyncDatalogr)

spdata <- Coords
spdataCoord <- spdata[order(spdata$region),]
disttab <- as.matrix(dist(spdataCoord[,2:3], upper=TRUE, diag=FALSE)/1000) # divide with 1000 to get km
colnames(disttab) <- spdataCoord$region
rownames(disttab) <- spdataCoord$region

SynchTab <- expand.grid(StartRho0=c(0.2,0.4,0.6,0.8), 
                        StartRhoInf=c(0,0.2,0.4,0.6), 
                        StartLambda=log(c(10,20,40,80,160,320, 1000)), 
                        species=Species, 
                        SynchVar="log.r", 
                        MaxLik=NA,
                        rho.0=NA, rho.inf=NA, Lambda=NA, MeanRhos=NA, N.areas=NA,
                        stringsAsFactors=FALSE)

SynchTab <- SynchTab[SynchTab$StartRho0>=SynchTab$StartRhoInf,] # starting values for rho0 always greater than rhoinf

for(i in 1:nrow(SynchTab))
{
  start <- SynchTab[i,c("StartRho0", "StartRhoInf", "StartLambda")]
  
  sp <- SynchTab[i,"species"]
  svar <- SynchTab[i,"SynchVar"]
  
  synchdatawide <- SyncLists[[svar]][[sp]]
  areas <- substr(colnames(synchdatawide), nchar(svar)+2, nchar(colnames(synchdatawide))) # need to get the names of areas with data, not all areas have all birds
  disttab2 <- disttab[areas, areas] #distance matrix must have same number of row and columns as number of areas
  gausssynch <- optim(start, GaussSyncF, distmat=disttab2, x=synchdatawide)
  gausssynchFix <- optimize(GaussSyncFFix, c(-1,1), distmat=disttab2, x=synchdatawide)
  SynchTab[i,"MaxLik"] <- gausssynch$value
  SynchTab[i,"rho.0"] <- exp(gausssynch$par[1])/(1+exp(gausssynch$par[1]))
  SynchTab[i,"rho.inf"] <- SynchTab[i,"rho.0"]*exp(gausssynch$par[2])/(1+exp(gausssynch$par[2]))
  SynchTab[i,"Lambda"] <- exp(gausssynch$par[3])
  SynchTab[i,"MeanRhos"] <- gausssynchFix$minimum
  SynchTab[i,"N.areas"] <- length(areas)
  print(paste(i, "of", nrow(SynchTab)))
  flush.console()
}

filename= paste0(country,"_", taxa, "_", method)
st=format(Sys.time(), "%d-%m-%Y")

synchtabb <- with(SynchTab, tapply(MaxLik, list(species, SynchVar), min))
SynchTabB <- SynchTab[SynchTab$MaxLik%in%synchtabb,]
SynchTabB<-distinct(SynchTabB, species, SynchVar, .keep_all=TRUE)
SynchTabB['Dataset'] = country
SynchTabB['Taxa'] = taxa


#### Calculate mean correlations and Fisher z-transform
SpeciesMeanbyDistance<-data.frame(species=character(0), Mean250kmZ=numeric(0), meandist=numeric(0), mediandist=numeric(0))
Species <- (unique(SynchTab$species))
for(i in 1:length(Species))
{
  sp <- Species[i]
  synchdatawider <- SyncLists[["log.r"]][[sp]]
  areas <- substr(colnames(synchdatawider), nchar(svar)+2, nchar(colnames(synchdatawider))) # need to get the names of areas with data, not all areas have all birds
  disttab2 <- disttab[areas, areas] #distance matrix must have same numbe
  TFmatrix <- disttab2<251
  popsunder250<-disttab2*TFmatrix
  popsunder250[popsunder250==0] <- NA
  meandist <- mean(popsunder250, na.rm=TRUE)
  mediandist <- median(popsunder250, na.rm=TRUE)
  cormatr <- rcorr(synchdatawider)
  cormatr$r[cormatr$r==1] <- NA
  FisherCorMat<-FisherZ(cormatr$r)
  FisherCorMat2<- lower.tri(FisherCorMat, diag = FALSE)
  FisherCorMat3 <- FisherCorMat*FisherCorMat2
  FisherCorMat3[FisherCorMat3==0] = NA
  FisherCorMat3[FisherCorMat3=="NaN"] = NA
  mean250r<-mean(FisherCorMat3[disttab2<251], na.rm=TRUE) 
  FisherZInvCorr<- FisherZInv(mean250r)
  TFTEST1<-lower.tri(cormatr$r)
  cormatlowertri<-cormatr$r*TFTEST1
  cormatlowertri[cormatlowertri==0] <- NA
  cormatdist <- cormatlowertri * TFmatrix
  cormatdist[cormatdist==0] <- NA
  mean250km <- mean(cormatlowertri[disttab2<251], na.rm=TRUE) 
  
  newrow<- data.frame(species = as.character(sp), Mean250kmZ = FisherZInvCorr, meandist = meandist, mediandist = as.numeric(mediandist), na.rm=TRUE) 
  SpeciesMeanbyDistance <- rbind(SpeciesMeanbyDistance, newrow)              # Apply rbind function
  
} 
meandistances <- SpeciesMeanbyDistance %>% dplyr::select("species","Mean250kmZ","mediandist", "meandist")
