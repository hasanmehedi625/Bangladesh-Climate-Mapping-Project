install.packages(c("raster", "terra", "ggplot2", "gstat", "fields", "stars"))
install.packages('ncdf4')
rm(list = ls())

library(raster)
library(terra)
library(ggplot2)
library(gstat)
library(fields)
library(stars)

setwd("E:\\Master's DSCR Resources\\507-Advanced RS and GIS\\Fahim bhai\\5efd9210f37768445cb780eeaccca9a0")
temperature_data <- rast("data_0.nc")

# Calculate the mean temperature across the time dimension
mean_temperature <- app(temperature_data, fun = mean)

# Convert to data frame for ggplot
temperature_df <- as.data.frame(mean_temperature, xy = TRUE)
temperature_data<- rast(temperature_df)

# Check the first few rows and column names of the data frame
head(temperature_df)
colnames(temperature_df)

# Plot the data using ggplot2
# Assuming the column with temperature values is named 'temperature'
ggplot(temperature_df) +
  geom_tile(aes(x = x, y = y, fill = mean)) +  # Use the correct column name here
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Mean 2m Temperature for April (2010-2024)",
       x = "Longitude", y = "Latitude", fill = "Temperature (°C)")

ggsave("mean_temperature_map_full_image.png", width = 8, height = 6)

# Read the NetCDF file
netcdf_data <- rast("data_0.nc")
#nc_data <- brick("data_0.nc")

# Convert the raster object to a data frame
df <- as.data.frame(netcdf_data, xy = TRUE)

# Check the structure of the data frame (optional)
head(df)

# change heading
colnames(df) <- c("Lat", 'Log', "2010_April", "2011_April", "2012_April", "2013_April", "2014_April", "2015_April", "2016_April", "2017_April", "2018_April", "2019_April", "2020_April", "2021_April", "2022_April", "2023_April", "2024_April")

head(df)
# Write the data frame to a CSV file
write.csv(df, "output_data_full_data.csv", row.names = FALSE)

# Load necessary libraries
library(terra)
library(sf)
#library(ggplot2)

# Set working directory
#setwd("J:\\Fahim Vai\\5efd9210f37768445cb780eeaccca9a0")

# Read the NetCDF file (temperature data in Kelvin)
temperature_data <- rast("data_0.nc")

# Calculate the mean temperature across the time dimension
mean_temperature_kelvin <- app(temperature_data, fun = mean)

# Convert the mean temperature from Kelvin to Celsius
mean_temperature_celsius <- mean_temperature_kelvin - 273.15

# Read the Bangladesh shapefile
bangladesh_shapefile <- st_read("BD_Zillas.shp")  # Replace with correct path

# Reproject the Bangladesh shapefile to match the raster CRS (if necessary)
bangladesh_shapefile <- st_transform(bangladesh_shapefile, crs(temperature_data))

# Convert the shapefile to a 'SpatVector' object for compatibility with 'terra'
bangladesh_sf <- vect(bangladesh_shapefile)

# Clip the raster (mean temperature in Celsius) using the Bangladesh shapefile
clipped_raster_celsius <- mask(mean_temperature_celsius, bangladesh_sf)

# Convert the clipped raster to a data frame for ggplot
clipped_df_celsius <- as.data.frame(clipped_raster_celsius, xy = TRUE)

# Check the first few rows of the data frame
head(clipped_df_celsius)

# Plot the clipped raster (mean temperature in Celsius) using ggplot2
ggplot(clipped_df_celsius) +
  geom_tile(aes(x = x, y = y, fill = mean)) +  # Replace 'layer' with the correct column name
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Temporal Heat Map of Bangladesh: Mean 2m Temperature Trends (2010–2024)",
       x = "Longitude", y = "Latitude", fill = "Temperature (°C)")

# Save the plot as a PNG file
ggsave("Temporal_Heat_Map_of_Bangladesh_Mean_2m_Temperature_Trends_(2010–2024).png", width = 4, height = 3, bg = "white")
theme(plot.title = element_text(size = 10)) 
# Optional: Save the clipped raster in Celsius as a GeoTIFF file
writeRaster(clipped_raster_celsius, "Mean_temperature_data_celsius_BD.tif", overwrite = TRUE)

df_bd <- as.data.frame(clipped_df_celsius, xy = TRUE)
#df_bd <- as.data.frame(clipped_raster_celsius, xy = TRUE)

# Check the structure of the data frame (optional)
head(df_bd)

# Write the data frame to a CSV file
write.csv(df_bd, "output_data_BD.csv", row.names = FALSE)

