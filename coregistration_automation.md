Drone Multispectral Coregistration & Visualization
================
Iris Nana Obeng
2026-02-04

``` r
library(terra)

setwd("~/Desktop/drone ms coregistration")
```

# 2. Load multispectral bands (Red = reference)

``` r
files <- c(
  "DJI_20251115110614_0003_MS_G (1).tiff",
  "DJI_20251115110614_0003_MS_R.tiff",
  "DJI_20251115110614_0003_MS_NIR.tiff"
)

ms <- rast(files)
```

    ## Warning: [rast] unknown extent

``` r
names(ms) <- c("Green", "Red", "NIR")
```

# 3. Alignment functions (matrix-based correlation)

``` r
align_band_simple <- function(ref_mat, target_mat, max_shift = 20) {
  best_cor <- -Inf; best_dx <- 0; best_dy <- 0
  rows <- nrow(ref_mat); cols <- ncol(ref_mat)
  
  for (dx in -max_shift:max_shift) {
    for (dy in -max_shift:max_shift) {
      shifted <- matrix(NA, rows, cols)
      r_dst <- (1:rows) + dy
      c_dst <- (1:cols) + dx
      ok_r <- r_dst >= 1 & r_dst <= rows
      ok_c <- c_dst >= 1 & c_dst <= cols
      shifted[r_dst[ok_r], c_dst[ok_c]] <-
        target_mat[(1:rows)[ok_r], (1:cols)[ok_c]]
      
      overlap <- !is.na(ref_mat) & !is.na(shifted)
      if (sum(overlap) > 1000) {
        cc <- cor(ref_mat[overlap], shifted[overlap], use = "complete.obs")
        if (!is.na(cc) && cc > best_cor) {
          best_cor <- cc; best_dx <- dx; best_dy <- dy
        }
      }
    }
  }
  list(dx = best_dx, dy = best_dy, cor = best_cor)
}
```

``` r
apply_shift <- function(mat, dx, dy) {
  rows <- nrow(mat); cols <- ncol(mat)
  shifted <- matrix(NA, rows, cols)
  r_dst <- (1:rows) + dy
  c_dst <- (1:cols) + dx
  ok_r <- r_dst >= 1 & r_dst <= rows
  ok_c <- c_dst >= 1 & c_dst <= cols
  shifted[r_dst[ok_r], c_dst[ok_c]] <-
    mat[(1:rows)[ok_r], (1:cols)[ok_c]]
  shifted
}
```

# 4. Align Green and NIR to Red

``` r
ref_red_mat <- as.matrix(ms$Red, wide = TRUE)
aligned_rasters <- list(Red = ms$Red)

for (b in c("Green", "NIR")) {
  target_mat <- as.matrix(ms[[b]], wide = TRUE)
  res <- align_band_simple(ref_red_mat, target_mat)
  shifted_mat <- apply_shift(target_mat, res$dx, res$dy)
  
  r_aligned <- rast(shifted_mat)
  ext(r_aligned) <- ext(ms$Red)
  crs(r_aligned) <- crs(ms$Red)
  
  aligned_rasters[[b]] <- r_aligned
  
  cat(sprintf("✓ %s aligned (dx=%d, dy=%d, cor=%.3f)\n",
              b, res$dx, res$dy, res$cor))
}
```

    ## ✓ Green aligned (dx=20, dy=15, cor=0.457)
    ## ✓ NIR aligned (dx=-16, dy=-20, cor=-0.024)

# 5. RGB stacks for plotting

``` r
# Red + Green
rg_before <- c(ms$Red, ms$Green, ms$Green)
rg_after  <- c(aligned_rasters$Red, aligned_rasters$Green, aligned_rasters$Green)

# Red + NIR
rn_before <- c(ms$Red, ms$NIR, ms$NIR)
rn_after  <- c(aligned_rasters$Red, aligned_rasters$NIR, aligned_rasters$NIR)

# Green + Red + NIR
rgb_before <- c(ms$Green, ms$Red, ms$NIR)
rgb_after  <- c(aligned_rasters$Green, aligned_rasters$Red, aligned_rasters$NIR)
```

# 6. Plotting

``` r
# Red + Green
par(mfrow = c(1,2), mar = c(2,2,3,1))
plotRGB(rg_before, r=1, g=2, b=3, stretch="lin",
        main="Red + Green\nBefore Alignment")
plotRGB(rg_after, r=1, g=2, b=3, stretch="lin",
        main="Red + Green\nAfter Alignment")
```

<img src="coregistration_automation_files/figure-gfm/plotting-1.png" style="display: block; margin: auto;" />

``` r
par(mfrow = c(1,1))

# Red + NIR
par(mfrow = c(1,2), mar = c(2,2,3,1))
plotRGB(rn_before, r=1, g=2, b=3, stretch="lin",
        main="Red + NIR\nBefore Alignment")
plotRGB(rn_after, r=1, g=2, b=3, stretch="lin",
        main="Red + NIR\nAfter Alignment")
```

<img src="coregistration_automation_files/figure-gfm/plotting-2.png" style="display: block; margin: auto;" />

``` r
par(mfrow = c(1,1))

# Green + Red + NIR
par(mfrow = c(1,2), mar = c(2,2,3,1))
plotRGB(rgb_before, r=3, g=2, b=1, stretch="lin",
        main="Green + Red + NIR\nBefore Alignment")
plotRGB(rgb_after, r=3, g=2, b=1, stretch="lin",
        main="Green + Red + NIR\nAfter Alignment")
```

<img src="coregistration_automation_files/figure-gfm/plotting-3.png" style="display: block; margin: auto;" />

``` r
par(mfrow = c(1,1))
```

``` r
cat("✓ Coregistration and all classic RGB combinations plotted successfully.\n")
```

    ## ✓ Coregistration and all classic RGB combinations plotted successfully.
