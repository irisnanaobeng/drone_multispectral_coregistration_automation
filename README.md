# Drone Multispectral Coregistration (Crop Dataset)

## Overview

This project implements a complete automatic band coregistration workflow for drone-based multispectral imagery using R.

The goal is to spatially align:

- Green band  
- Red band (reference)  
- Near-Infrared (NIR) band  

into a geometrically consistent stack suitable for band fusion and vegetation analysis.

The Red band is used as the spatial reference, while the Green and NIR bands are aligned to it.


## Methodology

### 1. Data Loading

Multispectral GeoTIFF files are loaded using the `terra` package:

- `ms_crop_gre.tif`
- `ms_crop_red.tif`
- `ms_crop_nir.tif`

These are stacked into a multi-layer raster object for processing.


### 2. Green Band Alignment (Integer Pixel Shift)

The Green band is aligned to the Red band using:

- Brute-force search over integer pixel shifts  
- Pearson correlation maximization  
- Overlapping valid pixel masking  

The best spatial shift `(dx, dy)` is then applied using matrix translation.


### 3. NIR Band Alignment (Edge-Based)

Due to spectral differences between Red and NIR, direct intensity matching is less effective.

Instead, alignment is performed using:

- Sobel edge detection  
- Edge-based correlation maximization  
- Integer pixel shift estimation  

This improves structural alignment robustness between bands.


### 4. Sub-Pixel Refinement

After coarse alignment, sub-pixel refinement is performed using:

- `imshift()` from the `imager` package  
- Local search around the optimal integer shift  
- Correlation maximization on overlapping pixels  

This achieves approximately 0.1-pixel alignment precision.


### 5. Visualization

Aligned stacks are visualized using `plotRGB()`:

- Red + Green  
- Red + NIR  
- Green + Red + NIR composite  

Before and after alignment comparisons are provided.


## Requirements

Required R packages:

"terra", "imager", "knitr"

## Output

The workflow produces:

Aligned Green raster

Aligned and refined NIR raster

RGB visual comparison panels

Correlation statistics for each alignment stage

## Author
Iris Nana Obeng
MSc. Global Change Ecology and Sustainable Development
