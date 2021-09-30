###################################################################################################
### sm stations to sm met.R
###
### Function to covert individual SnowModel station to a SnowModle Met File.
###
### Note:
###   - This code will make the standard met file, sw extra met file, and lw extra met file
###
###################################################################################################


#### Loading packages #############################################################################
#These must be already installed on your system 
library(dplyr)      # data manipulation
library(devtools)   # developer tools, simplifying tasks

### setting working file paths ####################################################################
sfl <- dirname(rstudioapi::getActiveDocumentContext()$path) # This is the filepath for this R script

### Assinging a data time range ###################################################################
start <- as.POSIXct("2015-10-01 00:00", tz = "UTC")
end <- as.POSIXct("2016-09-30 23:00", tz = "UTC")


### met files #####################################################################################

# Listing met stations files
station_files <- list.files(paste0(sfl,"/SnowModel Stations"),pattern=".csv", full.names=TRUE)

# counting number of met stations and timesteps
station_count <- length(station_files)
nt <- length(seq.POSIXt(from = start, to = end, by = "1 hour"))

# Loading met stations to a list
Stations <- list()
for (i in 1:length(station_files)){
  Stations[[i]] <-  read.csv(station_files[i])
  Stations[[i]]$datetime <- as.POSIXct(Stations[[i]]$datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")
}

# Filtering by station date range
for (i in 1:length(station_files)){
  Stations[[i]] <- Stations[[i]] %>% filter(datetime >= start, datetime <= end)
}

# Removing datetime element
for (i in 1:length(station_files)){
  Stations[[i]] <- Stations[[i]][2:14]
}

# Building the SM met file
# Note: this will take some time to run
data <- c(station_count,rep(NA, 12))     # start of the return data list
for (i in 1:nt){
  data <- rbind(data, do.call(rbind,(lapply(Stations, function(x) x[i,]))))
  data <- rbind(data,c(station_count,rep(NA, 12)))
}

# saving the sm file to the to the SnowModel_met folder
write.table(data,
            file = paste0(sfl,"/SnowModel Met/Sagehen_Met_WY2016_V1.dat"),
            append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

# clearing the data
rm(data,Stations,station_files,station_count,nt)



### SW files ######################################################################################
# Listing met stations files
station_files <- list.files(paste0(sfl,"/SnowModel Stations SW"),pattern=".csv", full.names=TRUE)

# counting number of met stations and timesteps
station_count <- length(station_files)
nt <- length(seq.POSIXt(from = start, to = end, by = "1 hour"))

# Loading met stations to a list
Stations <- list()
for (i in 1:length(station_files)){
  Stations[[i]] <-  read.csv(station_files[i])
  Stations[[i]]$datetime <- as.POSIXct(Stations[[i]]$datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")
}

# Filtering by station date range
for (i in 1:length(station_files)){
  Stations[[i]] <- Stations[[i]] %>% filter(datetime >= start, datetime <= end)
}

# Removing datetime element
for (i in 1:length(station_files)){
  Stations[[i]] <- Stations[[i]][2:10]
}

# Building the SM met file
# Note: this will take some time to run
data <- c(station_count,rep(NA, 12))     # start of the return data list
for (i in 1:nt){
  data <- rbind(data, do.call(rbind,(lapply(Stations, function(x) x[i,]))))
  data <- rbind(data,c(station_count,rep(NA, 12)))
}

# saving the sm file to the to the SnowModel_met folder
write.table(data,
            file = paste0(sfl,"/SnowModel Met SW/Sagehen_SW_WY2016_V1.dat"),
            append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

# clearing the data
rm(data,Stations,station_files,station_count,nt)


### LW files ######################################################################################
# Listing met stations files
station_files <- list.files(paste0(sfl,"/SnowModel Stations LW"),pattern=".csv", full.names=TRUE)

# counting number of met stations and timesteps
station_count <- length(station_files)
nt <- length(seq.POSIXt(from = start, to = end, by = "1 hour"))

# Loading met stations to a list
Stations <- list()
for (i in 1:length(station_files)){
  Stations[[i]] <-  read.csv(station_files[i])
  Stations[[i]]$datetime <- as.POSIXct(Stations[[i]]$datetime, format="%Y-%m-%d %H:%M:%S", tz="UTC")
}

# Filtering by station date range
for (i in 1:length(station_files)){
  Stations[[i]] <- Stations[[i]] %>% filter(datetime >= start, datetime <= end)
}

# Removing datetime element
for (i in 1:length(station_files)){
  Stations[[i]] <- Stations[[i]][2:10]
}

# Building the SM met file
# Note: this will take some time to run
data <- c(station_count,rep(NA, 12))     # start of the return data list
for (i in 1:nt){
  data <- rbind(data, do.call(rbind,(lapply(Stations, function(x) x[i,]))))
  data <- rbind(data,c(station_count,rep(NA, 12)))
}

# saving the sm file to the to the SnowModel_met folder
write.table(data,
            file = paste0(sfl,"/SnowModel Met LW/Sagehen_LW_WY2016_V1.dat"),
            append = FALSE, quote = TRUE, sep = "\t",
            eol = "\n", na = "", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

# clearing the data
rm(data,Stations,station_files,station_count,nt)

