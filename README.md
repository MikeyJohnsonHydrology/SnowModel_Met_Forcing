# SnowModel_Met_Forcing
Theses scritps can be used to read net CDF files from GRIDMET (Daily) or NLDAS2 (Hourly) and save the data into a format that can be read by SnowModel.

This process is run in two steps.

Step 1) Station information is saved independelty in the SnowModel format. Each station represnets a singel point in the gridded data set.
Step 2) All stations are combind into a singel SnowModel met file

Notes:
- Make sure to check file paths for loading and saving files
- Make sure to check for the spatial extents of your data before taking the time to read the files
- I typicaly save more data than I need for each station in Step 1 and then thin the data to the time I need in Step 2
- Depending on file size, you might need use alot of RAM for Step 2. Be aware of you computers RAM usage (I might try to rewrite this code at some point to fix  this problem).

GRIDMET Notes:
- All GRIDMET files are saved in a singel folder "grid_met_data", this included the GRIDMET elevation data

NDAS2 Notes:
- NLDAS2 files are available for each hour
- I have subset the data into folders year > month > day. The code below uses these file paths.
- The NLDAS2 data that I use is alos subset for California, you will neet to check the spatial location of the data you use.
