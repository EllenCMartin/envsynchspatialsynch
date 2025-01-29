################ CRU TS4 DATA GRID CLIMATE EXTRACTION #####
library(raster)
library(ncdf4)
lonlatproj <- "+proj=longlat +unit=dd +datum=WGS84"

## FILES NEEDED FOR THIS CODE TO RUN: 
## 1. "Coords": a long form file of X and Y coordinates of survey locations
##     "Coords" available for the following datasets and surveys used in analysis: 
## 1. National Finnish Butterfly Monitoring Scheme 
## 2. Sweden Songbirds, Finland Butterfly, UK Butterfly, Sweden Butterfly, Ireland Butterlly. All other countries survey location 

# Convert coordinates into the correct projection (lat-lon)
coordssf2 <- st_as_sf(as.data.frame(Coords), coords = c("X", "Y"), crs = 3035)
coordssf2 <- st_transform(coordssf2, crs = lonlatproj)
# Combine original coordinates with transformed ones
coords_df <- cbind(Coords, st_coordinates(coordssf2))
# Rename columns
colnames(coords_df) <- c("region", "Xutm", "Yutm", "X", "Y")
# Create a dataframe with only region, X, and Y
coords_df2 <- coords_df %>% dplyr::select(region, X, Y)
# Create a dataframe with X and Y coordinates only
coords_dataframe <- as.data.frame(cbind(coords_df$X, coords_df$Y))
# Use coords_dataframe (without region labels) or coords_df2 (with region labels) for CRU site extraction


setwd("CRUclimate")
nc.pre <- nc_open("2cru_ts4.06.1901.2021.pre.dat.nc")
nc.temp <- nc_open("2cru_ts4.06.1901.2021.tmp.dat.nc")
pre <- brick("2cru_ts4.06.1901.2021.pre.dat.nc", varname="pre")
temp <- brick("2cru_ts4.06.1901.2021.tmp.dat.nc", varname="tmp")
pre.sites <- data.frame(raster::extract(pre, coords_dataframe, ncol=2))
temp.sites <- data.frame(raster::extract(temp, coords_dataframe, ncol=2))

row.names(pre.sites) <- row.names(Coords)
row.names(temp.sites) <- row.names(Coords)
years <- 1901:2021
day <- c("01")
month <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
names(pre.sites) <- paste(rep(years, each=12), rep(month, times=121), rep(day, time=1452), sep="-")
names(temp.sites) <-paste(rep(years, each=12), rep(month, times=121), rep(day, time=1452), sep="-")

regions <- as.data.frame(coords_df2$region)
regions <- regions %>% 
  mutate(ID = 1:nrow(.))
names(regions)<- c("region", "ID")

## join these values with the site names 
pre.sites2 <- pre.sites %>% 
  mutate(ID = 1:nrow(.))

temp.sites2 <- temp.sites %>% 
  mutate(ID = 1:nrow(.))

Pre_with_site_labels <-  left_join(pre.sites2, regions)

Temp_with_site_labels <- left_join(temp.sites2, regions)

write.csv(Pre_with_site_labels, file="Precipitation Data.csv")
write.csv(Temp_with_site_labels, file="Temp Data.csv")

# Convert wide to long format for precipitation and temperature data
longprecip <- tidyr::gather(Pre_with_site_labels, date, precip, 1:1452)
longtemp <- tidyr::gather(Temp_with_site_labels, date, temp, 1:1452)

# Convert date format and extract month and year
longprecip <- longprecip %>% 
  dplyr::mutate(date = as.Date(date),
                month = format(date, "%m"),
                year = format(date, "%Y"))

longtemp <- longtemp %>% 
  dplyr::mutate(date = as.Date(date),
                month = format(date, "%m"),
                year = format(date, "%Y"))

# Join precipitation and temperature data together
environment <- longtemp %>% left_join(longprecip)

# Select data for summer months (April-November) from 1980 onwards
monthtarget <- c("04", "05", "06", "07", "08", "09", "10", "11")
environment5C_min <- filter(environment, year > 1980)
environment5C_minyears <- dplyr::filter(environment5C_min, month %in% monthtarget)

# Calculate monthly means for precipitation and temperature
monthlymeanprecip <- environment5C_minyears %>% 
  dplyr::group_by(region, year) %>% 
  dplyr::summarise(MeanPrecip = mean(precip))

monthlymeanprecipannual <- environment5C_min %>% 
  dplyr::group_by(region, year) %>% 
  dplyr::summarise(MeanPrecip = mean(precip))

monthlymeantemp <- environment5C_minyears %>% 
  dplyr::group_by(region, year) %>% 
  dplyr::summarise(MeanTemp = mean(temp))

monthlymeantempannual <- environment5C_min %>% 
  dplyr::group_by(region, year) %>% 
  dplyr::summarise(MeanTemp = mean(temp))

# Join monthly means for temperature and precipitation
monthlymeantempprecip <- as.data.frame(monthlymeanprecip %>% left_join(monthlymeantemp))
monthlymeantempprecipannual <- as.data.frame(monthlymeanprecipannual %>% left_join(monthlymeantempannual))

# Create wide-format data for precipitation and temperature
wideprecip <- monthlymeantempprecip %>% dplyr::select(region, MeanPrecip, year)
widetemp <- monthlymeantempprecip %>% dplyr::select(region, MeanTemp, year)
wideprecipannual <- monthlymeantempprecipannual %>% dplyr::select(region, MeanPrecip, year)
widetempannual <- monthlymeantempprecipannual %>% dplyr::select(region, MeanTemp, year)

# Convert to wide format using spread
wideprecip <- spread(wideprecip, key = year, value = MeanPrecip)
widetemp <- spread(widetemp, key = year, value = MeanTemp)
wideprecipannual <- spread(wideprecipannual, key = year, value = MeanPrecip)
widetempannual <- spread(widetempannual, key = year, value = MeanTemp)

# Ensure correct date format for plotting
environment5C_minyears$date <- as.Date(environment5C_minyears$date)
environment5C_minyears$month <- as.numeric(environment5C_minyears$month)

# Plot precipitation data by region
environment5C_minyears %>%  
  ggplot(aes(date, precip, col = region)) +
  geom_point(size = 1) +
  facet_wrap(~region, scales = 'free')

# Detrend precipitation data by region
for (i in 1:length(unique(environment5C_minyears$region))) {
  regionname <- unique(environment5C_minyears$region)[i]
  series <- filter(environment5C_minyears, region == regionname)
  detrendeddata <- detrend(series$precip, order = 1, lowess = FALSE, lowspan = 2/3)
  environment5C_minyears$detrended.precip[environment5C_minyears$region == regionname] <- detrendeddata
}

# Detrend temperature data by region
for (i in 1:length(unique(environment5C_minyears$region))) {
  regionname <- unique(environment5C_minyears$region)[i]
  series <- filter(environment5C_minyears, region == regionname)
  detrendeddata <- detrend(series$temp, order = 1, lowess = FALSE, lowspan = 2/3)
  environment5C_minyears$detrended.temp[environment5C_minyears$region == regionname] <- detrendeddata
}

# Calculate annual detrended means for temperature and precipitation
annualmeantemp.detrended <- environment5C_minyears %>% 
  dplyr::group_by(region, year) %>% 
  dplyr::summarise(MeanTemp.detrended = mean(detrended.temp))

annualmeanprecip.detrended <- environment5C_minyears %>% 
  dplyr::group_by(region, year) %>% 
  dplyr::summarise(MeanPrecip.detrended = mean(detrended.precip))

# Create wide-format data for detrended precipitation and temperature
widetempdetrended <- tidyr::spread(annualmeantemp.detrended, key = year, value = MeanTemp.detrended)
wideprecipdetrended <- tidyr::spread(annualmeanprecip.detrended, key = year, value = MeanPrecip.detrended)


