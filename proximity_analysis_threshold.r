# List of required packages
required_packages <- c("EBImage", "spatstat", "dplyr", "ggplot2")

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, dependencies=TRUE)

# Load libraries
library(EBImage)   
library(spatstat)  
library(dplyr)     
library(ggplot2)   

# Define image paths
image_folder <- readline(prompt = "Enter path to image folder: ")
output_folder <- readline(prompt = "Enter path for saving output images: ")

# Create output folder if it doesn't exist
if (!dir.exists(output_folder)) dir.create(output_folder)

# Generate file paths
n_images <- 24
control_images <- paste0(image_folder, control_prefix, 1:n_images, ".tif")
knockout_images <- paste0(image_folder, knockout_prefix, 1:n_images, ".tif")

# Combine into a single dataframe
image_data <- data.frame(
  genotype = c(rep("Control", 24), rep("Knockout", 24)),
  image_paths = c(control_images, knockout_images)
)

# Function to extract GFP and CD109 channels
extract_channels <- function(image_path) {
  img <- readImage(image_path)  # Load image
  
  if (length(dim(img)) == 3) {
    gfp_img <- channel(img, "green")  # GFP (green)
    cd109_img <- channel(img, "red")  # CD109 (red)
  } else {
    stop(paste("Image", image_path, "is not multichannel. Please check the format."))
  }
  
  return(list(gfp = gfp_img, cd109 = cd109_img, original = img))
}

# Function to analyze GFP/CD109 proximity with 22-pixel threshold
analyze_proximity <- function(image_path) {
  # Extract channels
  channels <- extract_channels(image_path)
  gfp_img <- channels$gfp
  cd109_img <- channels$cd109
  original_img <- channels$original

  # Binarize images for segmentation
  gfp_binary <- gfp_img > otsu(gfp_img)  
  cd109_binary <- cd109_img > otsu(cd109_img)

  # Label objects
  gfp_objects <- bwlabel(gfp_binary)
  cd109_objects <- bwlabel(cd109_binary)

  # Extract object centroids
  gfp_coords <- computeFeatures.moment(gfp_objects)[, c("m.cx", "m.cy"), drop=FALSE]
  cd109_coords <- computeFeatures.moment(cd109_objects)[, c("m.cx", "m.cy"), drop=FALSE]

  # If no objects detected, return NA
  if (nrow(gfp_coords) == 0 | nrow(cd109_coords) == 0) {
    return(data.frame(
      image = basename(image_path),
      total_GFP = nrow(gfp_coords),
      within_22px = NA,
      beyond_22px = NA
    ))
  }

  # Compute nearest neighbor distances
  gfp_points <- ppp(gfp_coords[,1], gfp_coords[,2], window=owin(c(0, dim(gfp_img)[1]), c(0, dim(gfp_img)[2])))
  cd109_points <- ppp(cd109_coords[,1], cd109_coords[,2], window=owin(c(0, dim(cd109_img)[1]), c(0, dim(cd109_img)[2])))
  
  distances <- nncross(gfp_points, cd109_points)$dist  # Nearest CD109 distance for each GFP object

  # Count objects based on 22-pixel threshold
  within_22px <- sum(distances <= 22, na.rm=TRUE)
  beyond_22px <- sum(distances > 22, na.rm=TRUE)

  # Create mask to overlay detected objects
  overlay_img <- original_img
  overlay_img[,,1] <- overlay_img[,,1] + (cd109_binary * 0.5)  # Red for CD109
  overlay_img[,,2] <- overlay_img[,,2] + (gfp_binary * 0.5)  # Green for GFP

  # Save overlay image
  output_path <- paste0(output_folder, basename(image_path))
  writeImage(overlay_img, output_path, quality = 85)

  return(data.frame(
    image = basename(image_path),
    total_GFP = nrow(gfp_coords),
    within_22px = within_22px,
    beyond_22px = beyond_22px
  ))
}

# Apply function to all images and compile results
proximity_results_threshold <- do.call(rbind, lapply(image_data$image_paths, analyze_proximity))

# Save results
results_file <- file.path(output_folder, "ProximityResultsThreshold.csv")
write.csv(image_data, file = results_file, row.names = FALSE)
# Print summary
print(proximity_results_threshold)
