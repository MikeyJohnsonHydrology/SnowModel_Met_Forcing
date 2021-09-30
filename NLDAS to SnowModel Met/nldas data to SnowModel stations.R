###################################################################################################
### nc data to sm stations.R
###
### Function to covert PRYSM corrects NLDAS data to SnowModel station format.
###
### Note:
###   - This data is provided by Sebastian Krogh
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

### setting working file paths ####################################################################
sfl <- dirname(rstudioapi::getActiveDocumentContext()$path) # This is the filepath for this R script

# Opening the dem.nc file
dem.nc <- nc_open(paste0(sfl,"/prism_corrected_nldas_Sierras/dem.nc"))

# Extracting spatial data, Each one of these point is a SnowModle station
grid_data <- expand.grid(x=1:32, y=1:48)
grid_data$lon <- dem.nc[["dim"]][["lon"]][["vals"]][grid_data$x]
grid_data$lat <- dem.nc[["dim"]][["lat"]][["vals"]][grid_data$y]
grid_data$elevation <- NA
grid_data$easting <- NA
grid_data$northing <- NA

for(i in 1: nrow(grid_data)){
  grid_data$elevation[i] <- ncvar_get(dem.nc,"Elev")[grid_data$x[i],grid_data$y[i]]
  grid_data$easting[i] <- as.numeric(round( lat_lon_to_UTM_10(grid_data$lat[i],grid_data$lon[i])[1] , digits = 4))
  grid_data$northing[i] <- as.numeric(round( lat_lon_to_UTM_10(grid_data$lat[i],grid_data$lon[i])[2] , digits = 4))
  }

grid_data <- grid_data %>% 
  filter(lat < 39.45) %>%     # Lat max
  filter(lat > 39.35) %>%     # Lat min
  filter(lon < -120.1) %>%   # Lon max
  filter(lon > -120.3)       # Lon min

grid_data <- grid_data[,c(1,2,5,6,7)] # Thinning the dataset


# Closing the dem.nc file
nc_close(dem.nc)
rm(dem.nc)


#Make data frame with datetime sequence for WY 2016
# Note: PST = Etc/GMT-8
time_data <- data.frame(datetime = seq.POSIXt(from = as.POSIXct("2015-10-01 00:00", tz='Etc/GMT+8'),
                                                to = as.POSIXct("2016-09-30 23:00", tz='Etc/GMT+8'),
                                                by = "1 hour"))

# Extract year, month, day and hour
time_data$year <- format(time_data$datetime, "%Y")
time_data$month <- format(time_data$datetime, "%m")
time_data$day <- format(time_data$datetime, "%d")
time_data$hour <- format(time_data$datetime, "%H")

# NLDAS uses the UTC (GTM) 8 hours difference from PST (Pacific Standard Time = Etc/GMT-8)
time_data$datetime_UTC <- format(time_data$datetime, tz="UTC")
time_data$year_UTC <- format(time_data$datetime, tz="UTC", "%Y")
time_data$month_UTC <- format(time_data$datetime, tz="UTC", "%m")
time_data$day_UTC <- format(time_data$datetime, tz="UTC", "%d")
time_data$hour_UTC <- format(time_data$datetime, tz="UTC", "%H")



### Loop throught stations for all time steps (Met Data)

#for(i in 1){ # Station Loop
for(i in 1:nrow(grid_data)){ # Station Loop
  
  tmp_station <- grid_data[i,] # Station information for the station loop
  tmp_sm_station <- NULL       # Station time series information that will be saved as a station file


  #for(j in 1:(365)){ # Timestep Loop
  for(j in 1:nrow(time_data)){ # Timestep Loop

    # Assigning a temporary time information file
    tmp_time <- time_data[j,]                       # Information for the timestep loop 
    tmp_j_loop_data <- cbind(tmp_station,tmp_time[,1:5])  # Station, time, and met data that will be added to tmp_sm_station

    # opening the temporary time step file
    tmp_nc <- nc_open(paste0(sfl,
                             "/prism_corrected_nldas_Sierras/",
                             tmp_time$year_UTC,
                             "/",
                             tmp_time$month_UTC,
                             "/",
                             tmp_time$day_UTC,
                             "/",
                             "prism_corrected_nldas_Sierras.",
                             tmp_time$year_UTC,
                             tmp_time$month_UTC,
                             tmp_time$day_UTC,
                             tmp_time$hour_UTC,
                             "00.nc"))

    # Tair (deg C)
    tmp_j_loop_data$Tair_K <- ncvar_get(tmp_nc,"TMP")[tmp_station$x,tmp_station$y]    # air temperature (K) ** at 2 meters above the surface

    # RH (%)
    tmp_j_loop_data$RH_degrees <- 0.263 * ncvar_get(tmp_nc,"PRES")[tmp_station$x,tmp_station$y] *  # surface pressure (Pa)
      ncvar_get(tmp_nc,"SPFH")[tmp_station$x,tmp_station$y] *                                      # specific humidity (kg/kg)
      (exp(17.67*(tmp_j_loop_data$Tair_K - 273.16)/(tmp_j_loop_data$Tair_K-29.65)))^-1             # eq: 0.263*p*q*(exp(17.67*(T-To)/(T-29.65)))^-1

    tmp_j_loop_data$RH_degrees <- ifelse(tmp_j_loop_data$RH_degrees > 100, 100, tmp_j_loop_data$RH_degrees) # Checking to make sure all values are less than 100%
    tmp_j_loop_data$RH_degrees <- ifelse(tmp_j_loop_data$RH_degrees < 0, 0, tmp_j_loop_data$RH_degrees)     # Checking to make sure all values are grater than 0%

    # Wind Speed (m/s)
    tmp_j_loop_data$wind_speed_mps <- sqrt(
      ncvar_get(tmp_nc,"UGRD")[tmp_station$x,tmp_station$y]^2 +   # U wind component (m/s) at 10 meters above the surface
        ncvar_get(tmp_nc,"VGRD")[tmp_station$x,tmp_station$y]^2   # V wind component (m/s) at 10 meters above the surface
    )

    # Wind dir (deg)
    tmp_j_loop_data$wind_dir_degrees <- acos(
      ncvar_get(tmp_nc,"VGRD")[tmp_station$x,tmp_station$y] / tmp_j_loop_data$wind_speed_mps
        ) * 180 / pi

    # Precip (mm/dt), note:(kg/m^2/dt) = (mm/dt)
    tmp_j_loop_data$precip_mmphr <- ncvar_get(tmp_nc,"APCP")[tmp_station$x,tmp_station$y]   # precipitation hourly total (kg/m^2) (mm/hr)


    # Addting the data to the dataframe to be saved
    tmp_sm_station <- rbind(tmp_sm_station,tmp_j_loop_data)

    # Closing and removing the current time step nc file
    nc_close(tmp_nc)
    rm(tmp_nc)

    # Removing the j-loop (timestep loop) specific tmp files
    rm(tmp_time, tmp_j_loop_data)

    }# End of the timestep (j-loop loop)

  # Creating the SM station dataframe to be saved
  tmp_sm_station <- tmp_sm_station %>% 
    mutate(Tair_C = Tair_K - 273.15,
           stn_id = i) %>%
    dplyr::select(datetime, year, month, day, hour, stn_id, easting, northing, 
                  elevation, Tair_C, RH_degrees, wind_speed_mps, wind_dir_degrees, 
                  precip_mmphr)

  # Saving the SM station dataframe
  write.table(tmp_sm_station,
              file = paste0(sfl,"/SM Met V1/SnowModel Stations/Station_",i,".csv"),
              sep=",",
              col.names=T,
              row.names=F)

  # Removing all the tmp data before the i-loop (station loop)
  rm(tmp_station, tmp_sm_station)

  }# End of the i-loop (station loop)

# SnowModel Format
# year   mo   dy    hr     stn_id  easting  northing  elevation   Tair     RH     speed    dir     precip
#(yyyy) (mm) (dd) (hh.hh) (number)   (m)       (m)      (m)        (C)    (%)     (m/s)   (deg)    (mm/dt)



### Loop throught stations for all time steps (SW & LW)

#for(i in 1){ # Station Loop
for(i in 1:nrow(grid_data)){ # Station Loop
  
  tmp_station <- grid_data[i,] # Station information for the station loop
  tmp_sm_station <- NULL       # Station time series information that will be saved as a station file
  

  for(j in 1:nrow(time_data)){ # Timestep Loop
    
    # Assigning a temporary time information file
    tmp_time <- time_data[j,]                       # Information for the timestep loop 
    tmp_j_loop_data <- cbind(tmp_station,tmp_time[,1:5])  # Station, time, and met data that will be added to tmp_sm_station
    
    # opening the temporary time step file
    tmp_nc <- nc_open(paste0(sfl,
                             "/prism_corrected_nldas_Sierras/",
                             tmp_time$year_UTC,
                             "/",
                             tmp_time$month_UTC,
                             "/",
                             tmp_time$day_UTC,
                             "/",
                             "prism_corrected_nldas_Sierras.",
                             tmp_time$year_UTC,
                             tmp_time$month_UTC,
                             tmp_time$day_UTC,
                             tmp_time$hour_UTC,
                             "00.nc"))
    
    # Shortwave Raidiation (Wpm2)
    tmp_j_loop_data$SWR_Wpm2 <- ncvar_get(tmp_nc,"DSWRF")[tmp_station$x,tmp_station$y] # surface downward shortwave radiation (W/m^2)
    
    # Longwave Raidiation (Wpm2)
    tmp_j_loop_data$LWR_Wpm2 <- ncvar_get(tmp_nc,"DLWRF")[tmp_station$x,tmp_station$y] # surface downward longwave radiation (W/m^2) 
    
    # Addting the data to the dataframe to be saved
    tmp_sm_station <- rbind(tmp_sm_station,tmp_j_loop_data)
    
    # Closing and removing the current time step nc file
    nc_close(tmp_nc)
    rm(tmp_nc)
    
    # Removing the j-loop (timestep loop) specific tmp files
    rm(tmp_time, tmp_j_loop_data)
    
  }# End of the timestep (j-loop loop)
  
  # Creating the SM station dataframe to be saved
  tmp_sm_station_sw <- tmp_sm_station %>% 
    mutate(stn_id = i) %>%
    dplyr::select(datetime, year, month, day, hour, stn_id, easting, northing, 
                  elevation, SWR_Wpm2)
  
  tmp_sm_station_lw <- tmp_sm_station %>% 
    mutate(stn_id = i) %>%
    dplyr::select(datetime, year, month, day, hour, stn_id, easting, northing, 
                  elevation, LWR_Wpm2)
  
  # Saving the SM station dataframe
  write.table(tmp_sm_station_sw,
              file = paste0(sfl,"/SM Met V1/SnowModel Stations SW/Station_SW_",i,".csv"),
              sep=",",
              col.names=T,
              row.names=F)
  
  write.table(tmp_sm_station_lw,
              file = paste0(sfl,"/SM Met V1/SnowModel Stations LW/Station_LW_",i,".csv"),
              sep=",",
              col.names=T,
              row.names=F)
  
  # Removing all the tmp data before the i-loop (station loop)
  rm(tmp_station, tmp_sm_station, tmp_sm_station_sw, tmp_sm_station_lw)
  
}# End of the i-loop (station loop)

# SW
# "yr" "mo" "dy" "hr" "stnID" "station_x" "station_y" "elevation" "incoming_SW_rad"

# LW
# "yr" "mo" "dy" "hr" "stnID" "station_x" "station_y" "elevation" "incoming_LW_rad"

