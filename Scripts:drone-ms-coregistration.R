---
  title: "DJI Multispectral Coregistration & Visualization"
author: "Iris Nana Obeng"
date: "2026-02-03"
output: github_document
---
  
  ```{r setup, include=FALSE}

# Setup: Loading required packages

library(terra)

# Set working directory
setwd("~/Desktop/drone ms coregistration")
```

## 1. Loading DJI multispectral bands

The dataset includes 3 bands: Green, Red, NIR.
The RED band will be used as reference for alignment.

```{r load-bands, message=FALSE}

# Load multispectral bands
files <- c(
  "DJI_20251115110614_0003_MS_G (1).tiff",
  "DJI_20251115110614_0003_MS_R.tiff",
  "DJI_20251115110614_0003_MS_NIR.tiff"
)
band_names <- c("Green", "Red", "NIR")

ms <- rast(files)

# Convert raster layers to matrices
band_matrices <- lapply(1:nlyr(ms), function(i) {
  as.matrix(ms[[i]], wide = TRUE)
})
names(band_matrices) <- band_names
```

## 2.Alignment function

Defining functions to shift bands relative to the reference RED band, maximizing correlation.

```{r align-bands}
align_band_simple <- function(ref_mat, target_mat, max_shift=20) {
  best_cor <- -Inf; best_dx <- 0; best_dy <- 0
  rows <- nrow(ref_mat); cols <- ncol(ref_mat)
  for(dx in -max_shift:max_shift) for(dy in -max_shift:max_shift) {
    shifted <- matrix(NA, rows, cols)
    r_dst <- (1:rows)+dy; c_dst <- (1:cols)+dx
    ok_r <- r_dst >= 1 & r_dst <= rows; ok_c <- c_dst >= 1 & c_dst <= cols
    shifted[r_dst[ok_r], c_dst[ok_c]] <- target_mat[(1:rows)[ok_r], (1:cols)[ok_c]]
    overlap <- !is.na(ref_mat) & !is.na(shifted)
    if(sum(overlap)>1000) {
      cc <- cor(ref_mat[overlap], shifted[overlap], use="complete.obs")
      if(!is.na(cc) && cc>best_cor) { best_cor <- cc; best_dx <- dx; best_dy <- dy }
    }
  }
  list(dx=best_dx, dy=best_dy, cor=best_cor)
}

apply_shift <- function(mat, dx, dy) {
  rows <- nrow(mat); cols <- ncol(mat); shifted <- matrix(NA, rows, cols)
  r_dst <- (1:rows)+dy; c_dst <- (1:cols)+dx
  ok_r <- r_dst >=1 & r_dst<=rows; ok_c <- c_dst>=1 & c_dst<=cols
  shifted[r_dst[ok_r], c_dst[ok_c]] <- mat[(1:rows)[ok_r], (1:cols)[ok_c]]
  shifted
}
```

## 2.1 Aligning Green and NIR bands to Red band (reference)

```{r align-green-nir}
ref_red <- band_matrices$Red
aligned_matrices <- list(Red=ref_red)
for(b in c("Green","NIR")) {
  res <- align_band_simple(ref_red, band_matrices[[b]])
  aligned_matrices[[b]] <- apply_shift(band_matrices[[b]], res$dx, res$dy)
}
```

## 2.2 Normalize and replace NA

Normalize values to [0,1] for plotting and replace NA values.

```{r normalize-replace-na} 
normalize <- function(x){ r<-range(x,na.rm=TRUE); (x-r[1])/(r[2]-r[1]) }
replace_na <- function(x){ x[is.na(x)]<-0; x }

R_b <- replace_na(normalize(band_matrices$Red))
G_b <- replace_na(normalize(band_matrices$Green))
N_b <- replace_na(normalize(band_matrices$NIR))

R_a <- replace_na(normalize(aligned_matrices$Red))
G_a <- replace_na(normalize(aligned_matrices$Green))
N_a <- replace_na(normalize(aligned_matrices$NIR))
```

## 3a. Plotting Helper Function (with Border & Caption)

```{r plot-function}
plot_rgb <- function(r, g, b, main_title, subtitle=NULL, border=TRUE) {
  arr <- array(c(r,g,b), dim=c(nrow(r),ncol(r),3))
  plot.new()
  rasterImage(arr,0,0,1,1)
  title(main_title, line=1.5)
  if(!is.null(subtitle)) mtext(subtitle, 3, 0.5, line=0.5)
  if(border) box(lwd=2)
}
```


## 3b.Visualization of original and aligned bands

### 3.1 Red + Green composite before and after alignment

```{r visualize-bands, fig.width=10, fig.height=5}

par(mfrow=c(1,2), mar=c(2,2,4,1))

par(mfrow=c(1,2), mar=c(2,2,3,1))
plot_rgb(R_b, G_b, 0, "Red + Green", "Before Alignment (Misaligned)")
plot_rgb(R_a, G_a, 0, "Red + Green", "After Alignment (Coregistered)")
par(mfrow=c(1,1))
```

### 3.2 Red + NIR composite before and after alignment

```{r visualize-red-nir, fig.width=10, fig.height=5}
par(mfrow=c(1,2), mar=c(2,2,3,1))
plot_rgb(R_b, 0, N_b, "Red + NIR", "Before Alignment (Misaligned)")
plot_rgb(R_a, 0, N_a, "Red + NIR", "After Alignment (Coregistered)")
par(mfrow=c(1,1))
```

### 3.3 Red + Green + NIR composite before and after alignment

```{r visualize-rgb, fig.width=10, fig.height=5}
par(mfrow=c(1,2), mar=c(2,2,3,1))
plot_rgb(N_b, R_b, G_b, "Red + Green + NIR", "Before Alignment (Misaligned)")
plot_rgb(N_a, R_a, G_a, "Red + Green + NIR", "After Alignment (Coregistered)")
par(mfrow=c(1,1))
```

## 4. High-Variance Region

Cropping a 300×300 patch with the highest variance to highlight misalignment more clearly.
```{r crop-high-variance}

high_variance_patch <- function(mat, size=300, step=50){
  rows<-nrow(mat); cols<-ncol(mat)
  best_var <- -Inf; best_rc <- c(1,1)
  for(r in seq(1, rows-size, by=step)) for(c in seq(1, cols-size, by=step)) {
    v <- var(as.vector(mat[r:(r+size-1), c:(c+size-1)]), na.rm=TRUE)
    if(!is.na(v) && v>best_var){ best_var<-v; best_rc<-c(r,c) }
  }
  best_rc
}

hv <- high_variance_patch(R_b)
r0<-hv[1]; c0<-hv[2]; crop<-function(x) x[r0:(r0+299), c0:(c0+299)]
```

### 4.1 Red + Green (High-Variance)
```{r visualize-high-variance, fig.width=10, fig.height=5}

par(mfrow=c(1,2), mar=c(2,2,3,1))
plot_rgb(crop(R_b), crop(G_b), 0, "Red + Green (High-Variance)", "Before Alignment (Misaligned)")
plot_rgb(crop(R_a), crop(G_a), 0, "Red + Green (High-Variance)", "After Alignment (Coregistered)")
par(mfrow=c(1,1))
```

### 4.2 Red + NIR (High-Variance)
```{r visualize-red-nir-high-variance, fig.width=10, fig.height=5}

par(mfrow=c(1,2), mar=c(2,2,3,1))
plot_rgb(crop(R_b), 0, crop(N_b), "Red + NIR (High-Variance)", "Before Alignment (Misaligned)")
plot_rgb(crop(R_a), 0, crop(N_a), "Red + NIR (High-Variance)", "After Alignment (Coregistered)")
par(mfrow=c(1,1))
```

### 4.3 Red + Green + NIR (High-Variance)

```{r visualize-rgb-high-variance, fig.width=10, fig.height=5}

par(mfrow=c(1,2), mar=c(2,2,3,1))
plot_rgb(crop(N_b), crop(R_b), crop(G_b), "Red + Green + NIR (High-Variance)", "Before Alignment (Misaligned)")
plot_rgb(crop(N_a), crop(R_a), crop(G_a), "Red + Green + NIR (High-Variance)", "After Alignment (Coregistered)")
par(mfrow=c(1,1))
```

```{r final-message, echo=FALSE}
cat("✓ All figures with borders and captions generated successfully.\n")
```
















