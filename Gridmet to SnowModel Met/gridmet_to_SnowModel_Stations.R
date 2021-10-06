###################################################################################################
### gridmet_SnowModel_Stations.R
###
### Function to covert GRIDMET data to SnowModel station format.
###
### By: Mikey Johnson
### Last Edited: 9/29/2021
###
###################################################################################################

#### Loading packages #############################################################################
#These must be already installed on your system 
library(dplyr)      # data manipulation
library(devtools)   # developer tools, simplifying tasks
library(ncdf4)      # netcdf manipulation
library(raster)     # rater manipulation
library(rgdal)      # geo-spatial data processing

#### Loading useful functions #####################################################################
source_url("https://raw.githubusercontent.com/MikeySnowHydro/Useful_R_Functions/master/lat_lon_to_UTM_10.R")
source_url("https://raw.githubusercontent.com/MikeySnowHydro/Useful_R_Functions/master/lat_lon_to_UTM_11.R")

### setting working file paths ####################################################################
sfl <- dirname(rstudioapi::getActiveDocumentContext()$path) # This is the filepath for this R script
setwd(sfl)


### Making a data frame for all the stations ######################################################

# Opening the dem.nc file
dem.nc <- nc_open(paste0(sfl,"/grid_met_data/metdata_elevationdata.nc"))
n_lat <- dim(ncvar_get(dem.nc,"elevation"))


# Extracting spatial data, Each one of these point is a SnowModle station
grid_data <- expand.grid(x=1:1386, y=1:585)
grid_data$lon <- dem.nc[["dim"]][["lon"]][["vals"]][grid_data$x]
grid_data$lat <- dem.nc[["dim"]][["lat"]][["vals"]][grid_data$y]
grid_data$elevation <- NA
grid_data$easting <- NA
grid_data$northing <- NA

# Filtering the gridmet data to a given the spatial extent
grid_data <- grid_data %>% 
  filter(lat < 37.20) %>%     # Lat max
  filter(lat > 37.10) %>%     # Lat min
  filter(lon < -118.90) %>%   # Lon max
  filter(lon > -119.10)       # Lon min


# Adding northing and easting data UTM Zone 11 
for(i in 1: nrow(grid_data)){
  grid_data$elevation[i] <- ncvar_get(dem.nc,"elevation")[grid_data$x[i],grid_data$y[i]]
  grid_data$easting[i] <- as.numeric(round( lat_lon_to_UTM_11(grid_data$lat[i],grid_data$lon[i])[1] , digits = 4))
  grid_data$northing[i] <- as.numeric(round( lat_lon_to_UTM_11(grid_data$lat[i],grid_data$lon[i])[2] , digits = 4))
}


# Thinning the data set (removing lat and lon)
grid_data <- grid_data %>% dplyr::select(x,y,elevation,easting,northing)


# Closing the dem.nc file
nc_close(dem.nc)
rm(dem.nc)


### Making a data frame for each calender year date #################################################

years <- 2019:2020 # (min year : max year)

# Make data frame with datetime sequence for each water year
# Note: PST = Etc/GMT-8
time_data <- data.frame(datetime = seq.POSIXt(from = as.POSIXct(paste0(min(years), "-01-01 00:00"), tz='Etc/GMT+8'),
                                              to = as.POSIXct(paste0(max(years),"-12-31 23:00"), tz='Etc/GMT+8'),
                                              by = "1 day"))
# Extract year, month, day and hour
time_data$year <- format(time_data$datetime, "%Y")
time_data$month <- format(time_data$datetime, "%m")
time_data$day <- format(time_data$datetime, "%d")
time_data$hour <- format(time_data$datetime, "%H")


### Loop throught stations for all time steps (Met Data) #############################################

#for(i in 1){ # Station Loop, i Loop
for(i in 1:nrow(grid_data)){ # Station Loop, i Loop
  
  tmp_station <- grid_data[i,] # Station information for the station loop
  tmp_sm_station <- NULL       # Station time series information that will be saved as a station file

  
  for(j in 1:length(years)){ # j loop, Timestep Loop
    tmp_yr <- years[j]
    tmp_ts_data <- time_data %>% filter(year == tmp_yr)
    tmp_ts_data$stn_id <- i
    tmp_ts_data$easting <- tmp_station$easting
    tmp_ts_data$northing <- tmp_station$northing
    tmp_ts_data$elevation <- tmp_station$elevation
    
    # Tair: Average of Tmax and Tmin (this code can be adjust this to improve results)
    Tmax.nc <- nc_open(paste0(sfl,"/grid_met_data/tmmx_",tmp_yr,".nc"))
    Tmin.nc <- nc_open(paste0(sfl,"/grid_met_data/tmmn_",tmp_yr,".nc"))
    
    tmp_ts_data$Tmax <- ncvar_get(Tmax.nc,"air_temperature")[tmp_station$x, tmp_station$y,]
    tmp_ts_data$Tmin <- ncvar_get(Tmin.nc,"air_temperature")[tmp_station$x, tmp_station$y,]
    
      # Calculating Tair
    tmp_ts_data <- tmp_ts_data %>% mutate(Tair=((Tmax+Tmin)/2)-273.15) # Mathematical Average & converting K to C
    tmp_ts_data$Tair <- round(tmp_ts_data$Tair, digits=2) 
    
      # Closing and removing the air temperature NC files
    nc_close(Tmax.nc)
    nc_close(Tmin.nc)
    rm(Tmax.nc,Tmin.nc)

        
    # Relative Humidity: Average of RHmax and RHmin (this code can be adjust this to improve results)
    RHmax.nc <- nc_open(paste0(sfl,"/grid_met_data/rmax_",tmp_yr,".nc"))
    RHmin.nc <- nc_open(paste0(sfl,"/grid_met_data/rmin_",tmp_yr,".nc"))
    
    tmp_ts_data$RHmax <- ncvar_get(RHmax.nc,"relative_humidity")[tmp_station$x, tmp_station$y,]
    tmp_ts_data$RHmin <- ncvar_get(RHmin.nc,"relative_humidity")[tmp_station$x, tmp_station$y,]

      # Calculating RH    
    tmp_ts_data <- tmp_ts_data %>% mutate(RH=(RHmax+RHmin)/2) # Mathematical Average
    
      # Closing and removing the relative humidity NC files
    nc_close(RHmax.nc)
    nc_close(RHmin.nc)
    rm(RHmax.nc,RHmin.nc)
  
    
    # Wind Speed and Direction
    speed.nc <- nc_open(paste0(sfl,"/grid_met_data/vs_",tmp_yr,".nc"))
    dir.nc <- nc_open(paste0(sfl,"/grid_met_data/th_",tmp_yr,".nc"))
    
    tmp_ts_data$speed <- ncvar_get(speed.nc,"wind_speed")[tmp_station$x, tmp_station$y,]
    tmp_ts_data$dir <- ncvar_get(dir.nc,"wind_from_direction")[tmp_station$x, tmp_station$y,]
    
      # Closing and removing the wind NC files
    nc_close(speed.nc)
    nc_close(dir.nc)
    rm(speed.nc,dir.nc)
    
    
    # Precipitation
    precip.nc <- nc_open(paste0(sfl,"/grid_met_data/pr_",tmp_yr,".nc"))

    tmp_ts_data$precip <- ncvar_get(precip.nc,"precipitation_amount")[tmp_station$x, tmp_station$y,]

      # Closing and removing the precip NC files
    nc_close(precip.nc)
    rm(precip.nc)

    
    # Cleaning the data for SnowModle station format  
    tmp_ts_data <- tmp_ts_data %>% dplyr::select(datetime,
                                                 year,
                                                 month,
                                                 day,
                                                 hour,
                                                 stn_id,
                                                 easting,
                                                 northing,
                                                 elevation,
                                                 Tair,
                                                 RH,
                                                 speed,
                                                 dir,
                                                 precip)

    # Saving the data to tmp_sm_stations  
    tmp_sm_station <- rbind(tmp_sm_station,tmp_ts_data)

    # Removing the j-loop (timestep loop) specific tmp files
    rm(tmp_ts_data)
    
  } # End of timestep Loop (J Loop)

  # Saving the SM station dataframe
    write.table(tmp_sm_station,
                file = paste0(sfl,"/SnowModel Station Data/Station_",i,".csv"),
                sep=",",
                col.names=T,
                row.names=F)
    
    # Removing all the tmp data before the i-loop (station loop)
    rm(tmp_station, tmp_sm_station)
    
}# End of the i-loop (station loop)
  
  
  
  
  
  