# Drone Multispectral Coregistration (Crop Dataset)

## Overview

This project implements an automatic band coregistration workflow for drone-based multispectral imagery using **R**.

The objective is to spatially align the following spectral bands:

* **Green**
* **Red** (reference band)
* **Near-Infrared (NIR)**

into a geometrically consistent stack suitable for band fusion.

The **Red band** is used as the spatial reference, while the **Green** and **NIR** bands are aligned to it. The workflow performs **coarse alignment followed by sub-pixel refinement**, and provides visual comparisons of the results.



## Data

The workflow uses three multispectral GeoTIFF files:

* `ms_crop_gre.tif`
* `ms_crop_red.tif`
* `ms_crop_nir.tif`

These files are loaded and stacked using the **terra** package.

The expected band order is:

```
1 = Green
2 = Red
3 = NIR
```



## Methodology

### 1. Data Loading

Multispectral TIFF files are loaded using the **terra** package and combined into a multi-layer raster object for processing.

```r
ms <- rast(files)
```

Band ordering is explicitly defined to ensure reproducibility across datasets.



### 2. Green Band Alignment (Intensity-Based)

The Green band is aligned to the Red reference band using an **intensity-based approach**:

* brute-force search across integer pixel shifts
* Pearson correlation maximization
* overlapping valid pixel masking

The optimal shift `(dx, dy)` is then applied using **bilinear interpolation**.

This step corrects the main spatial displacement between the two bands.



### 3. NIR Band Alignment (Edge-Based)

Direct intensity matching between Red and NIR bands is often unreliable due to strong spectral differences.

To address this, the workflow uses **edge-based registration**:

* Gaussian smoothing of the NIR band
* Sobel edge detection
* correlation maximization between edge images
* integer pixel shift estimation

Edge features provide more stable structural information for cross-spectral alignment.



### 4. Sub-Pixel Refinement

After coarse alignment, a **local sub-pixel search** refines the NIR alignment.

The refinement:

* searches around the coarse shift
* tests fractional pixel shifts
* maximizes correlation on overlapping pixels

This improves the registration accuracy beyond integer pixel resolution.



### 5. Visualization

Alignment quality is visually evaluated using `plotRGB()`.

The workflow generates the following comparisons:

**Red + Green**
Before and after alignment.

**Red + NIR**
Used to evaluate structural alignment between spectral bands.

**Green + Red + NIR composite**
Full multispectral stack before and after coregistration.

These visualizations help confirm the reduction of spatial ghosting between bands.



## Requirements

Required R package:

```
terra
```

The workflow runs using **base R and the terra package**.


## Output

The script produces:

* aligned **Green raster**
* aligned **NIR raster**
* RGB visual comparisons before and after alignment
* correlation statistics for each alignment stage

Example output includes:

```
Green aligned (dx, dy, correlation)
NIR coarse aligned (dx, dy, edge correlation)
NIR refined (dx, dy, correlation)
```



## Author

**Iris Nana Obeng**
MSc Global Change Ecology and Sustainable Development
