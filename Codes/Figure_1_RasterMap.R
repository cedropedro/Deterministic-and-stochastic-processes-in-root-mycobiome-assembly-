 # DESCRIPTION: Code for raster Map - FIGURE 1
# Modified from sergiocostafh: https://github.com/sergiocostafh/ggplot2_dataviz/blob/master/Raster%20maps%20with%20geom_raster().R


# Load required packages
library(raster)
library(ggplot2)
library(RColorBrewer)
library(sf)
library(maps)
library(mapdata)

# Retrieving data from worldclim - 12th layer is the precipitation, res. 0.5 requires lat/long
climate <- raster('wc2.1_30s_bio/wc2.1_30s_bio_12.tif')
coords_min <- extent(-98.5, -89.2, 37, 49.5)
climcrop <- crop(climate, coords_min)

# Define the plotting area limits
xlim = c(-98.5, -89.2)
ylim = c(37, 49.5)


# Custom color palette function
custom_palette <- colorRampPalette(c("#D22B2B","#D22B2B","red", "#FF8D21", "blue", "blue","darkblue","darkblue","navy"))
custom_palette <- colorRampPalette(c("#D22B2B","#D22B2B","red", "orange", "blue", "blue","darkblue","darkblue","navy"))

# Generate a palette with 15 colors
colors_15 <- custom_palette(40)

# Removing extra white space around the plot
par(mar = c(5, 4, 4, 5) + 0.1)

# Save the plot as a PNG or Pdf file without white space
#png("Figure_1_Sampling_Sites_Color.png", width = 5, height = 6, units = "in", res = 300) # pick this if you want a pdf
pdf("Figure_1_Sampling_Sites_Color.pdf", width = 5, height = 6)

# Create the main plot using plot and add legend title
plot(climcrop, col = colors_15, xlim = xlim, ylim = ylim,  
     xaxs = "i", yaxs = "i", 
     useRaster = TRUE, interpolate = FALSE, asp = NULL,
     xlab = expression(paste("Longitude (", degree, "W)")), 
     ylab = expression(paste("Latitude (", degree, "N)")),
     legend.args = list(text = "MAP (mm)", side = 4, line = 2.5, cex = 1))


# Define the states you want to overlay
states <- c("minnesota", "south dakota", "nebraska", "kansas")

# Overlay state borders for multiple states
for (state in states) {
  maps:: map("state", region = state, add = TRUE, lwd = 2, xlim = xlim, ylim = ylim)
}

# Sampling point
sampling_data <- read.csv('Site_coordinates.csv')
sampling_sf <- st_as_sf(sampling_data, coords = c("longitude", "latitude"), crs = 4326)

# Overlay sampling points
plot(st_geometry(sampling_sf), add = TRUE, col = "black", bg = "white", pch = 21, cex = 1.5)

# Define positions for each label (you can adjust these as needed)
positions <- c(4, 3, 3, 3, 2, 4, 3, 3, 3, 3) # Example positions for 10 points

# Adding labels for sampling points with non-overlapping settings
coords <- st_coordinates(sampling_sf)
for (i in 1:nrow(coords)) {
  text(coords[i, 1], coords[i, 2], labels = sampling_data$Code[i], cex = 0.7, font=2, pos = positions[i])
}

# Add state labels
state_labels <- data.frame(
  state = c("MN", "SD", "NE", "KS"),
  x = c(-94.9, -97.6, -97, -95.7),
  y = c(46, 45, 42, 39)
)

# Plot state labels
text(state_labels$x, state_labels$y, labels = state_labels$state, cex = 1.1, font = 2, col = "black")


dev.off()





