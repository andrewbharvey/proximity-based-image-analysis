# proximity-based-image-analysis

Proximity Analysis Script

Overview
This R script analyzes the spatial proximity between two cell types (GFP+ and CD109+) in multichannel microscopy images. It segments objects in the green (GFP) and red (CD109) channels, computes nearest-neighbor distances from GFP+ to CD109+ cells, and generates annotated overlay images and a summary of proximity results. This script should be modified for alternate applications when analyzing different fluorescent labels. The script also visualizes the proximity distributions between control and knockout genotypes using boxplots and jittered points.

The following R packages are automatically installed (if missing):

EBImage – for image reading and processing

spatstat – for spatial point pattern analysis

dplyr – for data manipulation

ggplot2 – for data visualization

How to Use
1. Set Up Your Images:
   
Prepare 24 images each for the "Control" and "Knockout" conditions. All images should be:

In TIFF format (.tif)

Multichannel, with GFP (or fluorophore of your choice) in the green channel and CD109 (or fluorophore of your choice) in the red channel

Named using a consistent pattern that includes the genotype prefix (e.g., Control1.tif, Knockout1.tif, etc.)

Set the appropriate values for control_prefix and knockout_prefix within the script. Example:

r
Copy
Edit
control_prefix <- "Control"
knockout_prefix <- "Knockout"

2. Run the Script:
   
When prompted:

Enter the full path to the image folder

Enter the desired output folder for processed images and results

Edit script to name fluorophores of choosing

The script will:

Load and segment each image

Compute the mean distance from each GFP+ object to the nearest CD109+ object

Save annotated overlays with labeled GFP/CD109 regions

Output a CSV file (ProximityResults.csv) with genotype and proximity values

Display a boxplot comparing the proximity across genotypes

Output:

Overlay Images: Saved in the output folder with GFP (green) and CD109 (red) masks.

Results File: ProximityResults.csv with the following columns:

genotype: Control or Knockout

image_paths: Path to the original image

proximity: Mean distance (in pixels) from GFP+ objects to their nearest CD109+ object for a single image

Boxplot: Visualization of proximity distributions between genotypes

Notes:

Images that do not have both GFP+ and CD109+ signals will be marked with NA for proximity.

Ensure all images are correctly formatted as multichannel TIFFs to avoid processing errors.
