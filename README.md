# VegChangeR

A comprehensive R package for vegetation change analysis using satellite imagery and raster data. Supports multiple time scales and memory-optimized processing for large datasets.

## Features

- **Multi-scale analysis**: 2-week, 1-month, 1-year, and 5-year change detection
- **Memory optimized**: Multiple processing modes for different RAM availability
- **Disk-based processing**: Handle very large datasets with chunk-wise processing
- **Seasonal consideration**: Smart date matching for period-over-period comparisons
- **Polygon analysis**: Extract and summarize changes for spatial features
- **Comprehensive visualization**: Built-in plotting and statistical summaries

## Installation

### From GitHub
```r
# Install devtools if not already installed
# install.packages("devtools")

devtools::install_github("DijoG/VegChangeR")
```
## Usage

### Quick start

```r
require(VegChangeR);require(terra)

# Load vegetation data (binary: 0 = non-vegetation, 1 = vegetation)
veg_data <- rast("path/to/vegetation_data.tif")
names(veg_data)

# Get memory recommendations
VegChangeR::recommend_memory_mode()
# "lowmem": 30% RAM - For systems with limited memory
# "medium": 50% RAM - Balanced performance
# "high": 70% RAM - Better performance
# "aggressive": 90% RAM - Maximum speed
# "auto": Let terra manage memory (default)

VegChangeR::monitor_memory()

# Calculate vegetation changes with stable category
changes <- VegChangeR::get_VC_optimized0(veg_data, "2023-06-15", processing_option = "medium")

# Clear memory cache
VegChangeR::clear_terra_cache()

# View results
plot(changes$oneM)  # 1-month changes

# Save results
VegChangeR::save_changes(changes, "path/to/output_directory")
```

### Memory-optimized processing

```r
# Load large vegetation data (binary: 0 = non-vegetation, 1 = vegetation)
large_veg_data <- rast("path/to/vegetation_data.tif")
names(veg_data)

# Run the memory-optimized process 
VegChangeR::CHUNKWISE_optimal_memfrac()
tictoc::tic()
VegChangeR::CHUNKWISE_memo_monitor()
changes <- VegChangeR::CHUNKWISE_get_VC_TOdisk(large_veg_data, "2025-06-15", temp_dir = "D:/temp_processing", auto_optimize = TRUE)
VegChangeR::CHUNKWISE_memo_monitor()
tictoc::toc()
```

### Polygon analysis
```r
# Load study area polygons
study_area <- st_read("path/to/polygons.shp")

# Extract changes to polygons
polygon_results <- VegChangeR::extract_changes_exact(changes, study_area)

# Analyze and visualize results
VegChangeR::analyze_vegetation_changes(polygon_results)
VegChangeR::plot_polygon_changes(polygon_results)
```





