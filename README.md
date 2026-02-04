# Drone Multispectral Coregistration & Visualization

Author: Iris Nana Obeng
Date: 2026-02-03


This repository contains an R script for coregistering DJI multispectral images (Green, Red, NIR) using a shift-based alignment method, and visualizing the results for both full-frame images and high-variance regions.

## Repository Overview

The project demonstrates:

 - Loading DJI multispectral bands into R (terra package).

 - Aligning Green and NIR bands to the Red band (used as reference)       using correlation-maximizing shifts.

 - Normalizing raster values and replacing missing pixels (NA).

 - Visualizing before and after coregistration for:

      Red + Green
      Red + NIR
      Red + Green + NIR
      High-variance regions (300×300 pixels)

Providing clear comparisons to highlight misalignments and improvements.

drone-ms-coregistration
  - ├─ DJI_20251115110614_0003_MS_G (1).tiff
  - ├─ DJI_20251115110614_0003_MS_R.tiff
  - ├─ DJI_20251115110614_0003_MS_NIR.tiff
  - ├─ DJI_3Band_Aligned.tif
  - ├─ Scripts/
  - │   └─ drone-ms-coregistration.R
  - ├─ Output/
  - └─ README.md

  - Scripts/ – Contains the R script for coregistration and visualization.

  

## Visualization 

Full-frame composites:

Red + Green
Red + NIR
Red + Green + NIR

High-variance region composites (to highlight misalignment):

Red + Green
Red + NIR
Red + Green + NIR

Each visualization shows before (misaligned) vs after (coregistered) alignment.

- Raw TIFF files – Input multispectral bands (Green, Red, NIR).

## Author
Iris Nana Obeng
(Msc. Global Change Ecology and Sustainable Development)
