#' Vegetation Change Analysis Tool
#' 
#' A comprehensive package for analyzing vegetation changes over multiple time scales
#' using satellite imagery and raster data. Supports memory-optimized processing
#' for large datasets.
#' 
#' @name VegChangeR
#' @docType package
#' @import terra sf dplyr tidyr exactextractr ggplot2 patchwork lubridate
NULL

#' Calculate vegetation changes across multiple time scales
#'
#' This function analyzes vegetation changes between a target date and multiple
#' comparison periods (2 weeks, 1 month, 1 year, 5 years) with memory-optimized
#' processing options.
#'
#' @param inputRAST A SpatRaster object containing vegetation data with dates in layer names
#' @param targetDATE Target date for analysis obtained from names of input raster stack
#' @param processing_option Memory processing option: "auto", "lowmem", "medium", 
#'        "high", or "aggressive". Default: "auto"
#'
#' @return A list containing:
#' \itemize{
#'   \item{twoW, oneM, oneY, fiveY: SpatRaster objects with vegetation changes}
#'   \item{metadata: List with processing information and date comparisons}
#' }
#'
#' @details
#' This function calculates vegetation changes where:
#' \itemize{
#'   \item{-1: Vegetation loss}
#'   \item{0: Stable vegetation} 
#'   \item{1: Vegetation gain}
#'   \item{NA: Non-vegetated areas or no data}
#' }
#'
#' The function includes seasonal consideration for year-based comparisons
#' and handles duplicate dates automatically.
#'
#' @examples
#' \dontrun{
#' # Load your vegetation raster
#' veg_data <- rast("path/to/vegetation_data.tif")
#' 
#' # Calculate changes for specific date
#' changes <- get_VC_optimized0(veg_data, "2023-06-15", processing_option = "medium")
#' 
#' # Plot results
#' plot(changes$oneY)
#' }
#'
#' @export
#' @import terra
get_VC_optimized0 <- function(inputRAST, 
                              targetDATE, 
                              processing_option = "auto") {
  
  # Set memory options based on processing choice
  if(processing_option == "medium") {
    terraOptions(memfrac = 0.5)
    cat("Low memory mode: Using 50% of available RAM\n")
  } else if(processing_option == "high") {
    terraOptions(memfrac = 0.7)
    cat("High processing mode: Using 70% of available RAM\n")
  } else if(processing_option == "lowmem") {
    terraOptions(memfrac = 0.3, tempdir = "D:/temp")
    cat("Very low memory mode: Using 30% of available RAM\n")
  } else if(processing_option == "aggressive") {
    terraOptions(memfrac = 0.9)
    cat("Aggressive mode: Using 90% of available RAM for maximum performance\n")
  } else {
    cat("Auto memory mode: Letting terra manage memory\n")
  }
  
  # Get the dates from layer names and handle duplicates
  layer_names = names(inputRAST)
  dates = as.Date(layer_names)
  unique_indices = !duplicated(dates)
  unique_rast = inputRAST[[unique_indices]]
  unique_dates = dates[unique_indices]
  unique_names = layer_names[unique_indices]
  
  # Sort by date
  date_order = order(unique_dates)
  sorted_rast = unique_rast[[date_order]]
  sorted_dates = unique_dates[date_order]
  sorted_names = unique_names[date_order]
  
  # Find the target date
  target_date = as.Date(targetDATE)
  target_idx = which(sorted_dates == target_date)
  
  if(length(target_idx) == 0) {
    stop("Target date ", targetDATE, " not found in raster layers")
  }
  
  cat("Target layer:", sorted_names[target_idx], "\n")
  
  # Function to find the closest available date with seasonal consideration
  find_seasonal_comparison = function(base_date, years_back = 0, approximate_days = NULL) {
    if (!is.null(approximate_days)) {
      target_lookup_date = base_date - approximate_days
    } else {
      target_lookup_date = base_date - years(years_back)
    }
    
    date_diffs = abs(sorted_dates - target_lookup_date)
    
    if (years_back > 0) {
      seasonal_window = 30
      in_season = abs(as.numeric(sorted_dates - target_lookup_date)) <= seasonal_window
      if (any(in_season)) {
        date_diffs[!in_season] = Inf
      }
    }
    
    closest_idx = which.min(date_diffs)
    
    return(list(
      idx = closest_idx, 
      date = sorted_dates[closest_idx],
      diff_days = as.numeric(abs(sorted_dates[closest_idx] - target_lookup_date)),
      name = sorted_names[closest_idx]
    ))
  }
  
  # Find comparison layers
  cat("Finding comparison dates...\n")
  twoW_info = find_seasonal_comparison(target_date, approximate_days = 14)
  oneM_info = find_seasonal_comparison(target_date, approximate_days = 30)
  oneY_info = find_seasonal_comparison(target_date, years_back = 1)
  fiveY_info = find_seasonal_comparison(target_date, years_back = 5)
  
  cat("Comparison layers found:\n")
  cat("2 weeks before:", twoW_info$name, "(difference:", twoW_info$diff_days, "days)\n")
  cat("1 month before:", oneM_info$name, "(difference:", oneM_info$diff_days, "days)\n")
  cat("1 year before:", oneY_info$name, "(difference:", oneY_info$diff_days, "days)\n")
  cat("5 years before:", fiveY_info$name, "(difference:", fiveY_info$diff_days, "days)\n")
  
  # Get the raster layers
  cat("Loading target and comparison layers...\n")
  target_layer = sorted_rast[[target_idx]]
  twoW_layer = sorted_rast[[twoW_info$idx]]
  oneM_layer = sorted_rast[[oneM_info$idx]]
  oneY_layer = sorted_rast[[oneY_info$idx]]
  fiveY_layer = sorted_rast[[fiveY_info$idx]]
  
  # Enhanced function to calculate vegetation change with stable vegetation
  calculate_change_with_stable = function(current, previous, period_name) {
    cat("Calculating", period_name, "changes with stable vegetation...\n")
    
    # Process in chunks for very large rasters (only in "high" and "lowmem" modes)
    if(processing_option %in% c("high", "lowmem") && ncell(current) > 1e7) {
      cat("  Large raster detected, processing in chunks...\n")
      
      change_rast = rast(current)
      min_blocks = ifelse(processing_option == "lowmem", 20, 10)
      bs = blockSize(change_rast, minblocks = min_blocks)
      
      for (i in 1:bs$n) {
        cat("  Processing chunk", i, "of", bs$n, "\n")
        
        # Read chunks from both rasters
        v_current = readValues(current, bs$row[i], bs$nrows[i])
        v_previous = readValues(previous, bs$row[i], bs$nrows[i])
        
        # STEP 1: Calculate gain/loss (-1, 1) through subtraction
        v_change = v_current - v_previous  # This gives -1, 0, 1
        
        # STEP 2: Calculate stable vegetation (2) through addition
        v_sum = v_current + v_previous  # This gives 0, 1, 2
        v_stable = ifelse(v_sum == 2, 2, 0)  # Keep only stable vegetation as 2
        
        # STEP 3: Combine both layers
        v_combined = v_change + v_stable  # This gives -1, 0, 1, 2, 3
        
        # STEP 4: Apply masking and final classification
        valid_mask = !is.na(v_current) & !is.na(v_previous)
        vegetation_mask = (v_current == 1) | (v_previous == 1)
        final_mask = valid_mask & vegetation_mask
        
        # Final classification:
        # -1: Loss (vegetation disappeared)
        # 0: No change (non-vegetated or transient)
        # 1: Gain (vegetation appeared)  
        # 2: Stable (vegetation persisted)
        v_final = v_combined
        v_final[!final_mask] = NA  # Mask out non-vegetated areas
        v_final[v_final == 0] = NA  # Mask out 0 values as requested
        v_final[v_final == 2] = 0   # Convert stable vegetation from 2 to 0
        v_final[v_final == 3] = 1   # Handle edge case where gain + stable = 3
        
        # Write chunk back to raster
        change_rast = writeValues(change_rast, v_final, bs$row[i])
      }
      
      change_rast = writeStop(change_rast)
      
    } else {
      # Standard processing for smaller rasters
      
      # STEP 1: Calculate gain/loss (-1, 1)
      change = current - previous
      
      # STEP 2: Calculate stable vegetation (2)  
      sum_layers = current + previous
      stable = classify(sum_layers, matrix(c(2, 2), ncol = 2))  # Keep only value 2 as 2
      stable = stable * 0 + ifel(sum_layers == 2, 2, 0)  # Alternative approach
      
      # STEP 3: Combine both layers
      combined = change + stable
      
      # STEP 4: Apply masking and final classification
      valid_mask = !is.na(current) & !is.na(previous)
      vegetation_mask = (current == 1) | (previous == 1)
      final_mask = valid_mask & vegetation_mask
      
      # Final classification
      change_final = combined * final_mask
      change_final[!final_mask] = NA  # Mask out non-vegetated areas
      change_final = ifel(change_final == 0, NA, change_final)  # Mask out 0 values
      change_final = ifel(change_final == 2, 0, change_final)   # Stable → 0
      change_final = ifel(change_final == 3, 1, change_final)   # Edge case
      
      change_rast = change_final
    }
    
    return(change_rast)
  }
  
  # Calculate changes for each time period
  cat("Calculating vegetation changes with stable category...\n")
  twoW_change = calculate_change_with_stable(target_layer, twoW_layer, "2-week")
  oneM_change = calculate_change_with_stable(target_layer, oneM_layer, "1-month")
  oneY_change = calculate_change_with_stable(target_layer, oneY_layer, "1-year")
  fiveY_change = calculate_change_with_stable(target_layer, fiveY_layer, "5-year")
  
  # Set names for the output rasters
  names(twoW_change) = paste0("twoW_change_", targetDATE)
  names(oneM_change) = paste0("oneM_change_", targetDATE)
  names(oneY_change) = paste0("oneY_change_", targetDATE)
  names(fiveY_change) = paste0("fiveY_change_", targetDATE)
  
  # Clear intermediate objects
  rm(target_layer, twoW_layer, oneM_layer, oneY_layer, fiveY_layer)
  gc()
  
  # Create output list
  result_list = list(
    twoW = twoW_change,
    oneM = oneM_change,
    oneY = oneY_change,
    fiveY = fiveY_change,
    metadata = list(
      target_date = target_date,
      target_name = sorted_names[target_idx],
      comparison_dates = list(
        twoW = twoW_info$date,
        oneM = oneM_info$date,
        oneY = oneY_info$date,
        fiveY = fiveY_info$date
      ),
      value_meaning = c(
        "-1" = "Vegetation loss",
        "0" = "Stable vegetation", 
        "1" = "Vegetation gain",
        "NA" = "Non-vegetated or no data"
      ),
      processing_mode = processing_option
    )
  )
  
  cat("Vegetation change calculation completed successfully!\n")
  cat("Final values: -1 (loss), 0 (stable), 1 (gain), NA (non-vegetated)\n")
  return(result_list)
}

#' Calculate vegetation changes (gain/loss only)
#'
#' Simplified version that only returns vegetation gain (+1) and loss (-1),
#' masking out stable vegetation and no-change areas.
#'
#' @inheritParams get_VC_optimized0
#'
#' @return A list containing SpatRaster objects with vegetation changes where:
#' \itemize{
#'   \item{-1: Vegetation loss}
#'   \item{1: Vegetation gain}
#'   \item{NA: Stable vegetation or non-vegetated areas}
#' }
#'
#' @examples
#' \dontrun{
#' changes <- get_VC_optimized(veg_data, "2023-06-15")
#' # Only shows gain (+1) and loss (-1)
#' }
#'
#' @export
get_VC_optimized <- function(inputRAST, 
                             targetDATE, 
                             processing_option = "auto") {
  
  # Set memory options based on processing choice - using your suggested fractions
  if(processing_option == "medium") {
    terraOptions(memfrac = 0.5)  # Use 50% of RAM
    cat("Low memory mode: Using 50% of available RAM\n")
  } else if(processing_option == "high") {
    terraOptions(memfrac = 0.7)  # Use 70% of RAM
    cat("High processing mode: Using 70% of available RAM\n")
  } else if(processing_option == "lowmem") {
    terraOptions(memfrac = 0.3, tempdir = "D:/temp")  # Use 30% RAM and custom temp dir
    cat("Very low memory mode: Using 30% of available RAM\n")
  } else if(processing_option == "aggressive") {
    terraOptions(memfrac = 0.9)  # Use 90% of RAM for maximum speed
    cat("Aggressive mode: Using 90% of available RAM for maximum performance\n")
  } else {
    # auto mode - let terra handle it (usually around 60-70%)
    cat("Auto memory mode: Letting terra manage memory\n")
  }
  
  # Get the dates from layer names and handle duplicates
  layer_names = names(inputRAST)
  dates = as.Date(layer_names)
  
  # Remove duplicate dates by keeping the first occurrence
  unique_indices = !duplicated(dates)
  unique_rast = inputRAST[[unique_indices]]
  unique_dates = dates[unique_indices]
  unique_names = layer_names[unique_indices]
  
  # Sort by date to ensure chronological order
  date_order = order(unique_dates)
  sorted_rast = unique_rast[[date_order]]
  sorted_dates = unique_dates[date_order]
  sorted_names = unique_names[date_order]
  
  # Find the target date
  target_date = as.Date(targetDATE)
  target_idx = which(sorted_dates == target_date)
  
  if(length(target_idx) == 0) {
    stop("Target date ", targetDATE, " not found in raster layers")
  }
  
  cat("Target layer:", sorted_names[target_idx], "\n")
  
  # Function to find the closest available date with seasonal consideration
  find_seasonal_comparison = function(base_date, years_back = 0, approximate_days = NULL) {
    if (!is.null(approximate_days)) {
      # For short-term comparisons (weeks, months)
      target_lookup_date = base_date - approximate_days
    } else {
      # For year-based comparisons - try to match the same seasonal period
      target_lookup_date = base_date - years(years_back)
    }
    
    # Find the closest date to the lookup date within a reasonable window
    date_diffs = abs(sorted_dates - target_lookup_date)
    
    # For year comparisons, prioritize dates within ±30 days of the target seasonal period
    if (years_back > 0) {
      seasonal_window = 30  # days
      in_season = abs(as.numeric(sorted_dates - target_lookup_date)) <= seasonal_window
      if (any(in_season)) {
        date_diffs[!in_season] = Inf  # Penalize dates outside seasonal window
      }
    }
    
    closest_idx = which.min(date_diffs)
    
    return(list(
      idx = closest_idx, 
      date = sorted_dates[closest_idx],
      diff_days = as.numeric(abs(sorted_dates[closest_idx] - target_lookup_date)),
      name = sorted_names[closest_idx]
    ))
  }
  
  # Find comparison layers with seasonal consideration
  cat("Finding comparison dates...\n")
  twoW_info = find_seasonal_comparison(target_date, approximate_days = 14)
  oneM_info = find_seasonal_comparison(target_date, approximate_days = 30)
  oneY_info = find_seasonal_comparison(target_date, years_back = 1)
  fiveY_info = find_seasonal_comparison(target_date, years_back = 5)
  
  cat("Comparison layers found:\n")
  cat("2 weeks before:", twoW_info$name, "(difference:", twoW_info$diff_days, "days)\n")
  cat("1 month before:", oneM_info$name, "(difference:", oneM_info$diff_days, "days)\n")
  cat("1 year before:", oneY_info$name, "(difference:", oneY_info$diff_days, "days)\n")
  cat("5 years before:", fiveY_info$name, "(difference:", fiveY_info$diff_days, "days)\n")
  
  # Get the raster layers
  cat("Loading target and comparison layers...\n")
  target_layer = sorted_rast[[target_idx]]
  twoW_layer = sorted_rast[[twoW_info$idx]]
  oneM_layer = sorted_rast[[oneM_info$idx]]
  oneY_layer = sorted_rast[[oneY_info$idx]]
  fiveY_layer = sorted_rast[[fiveY_info$idx]]
  
  # Memory-efficient function to calculate vegetation change
  calculate_change_efficient = function(current, previous, period_name) {
    cat("Calculating", period_name, "changes...\n")
    
    # Process in chunks for very large rasters (only in "high" and "lowmem" modes)
    if(processing_option %in% c("high", "lowmem") && ncell(current) > 1e7) {
      cat("  Large raster detected, processing in chunks...\n")
      
      # Calculate difference in chunks
      change_rast = rast(current)
      # Adjust chunks based on memory mode
      min_blocks = ifelse(processing_option == "lowmem", 20, 10)
      bs = blockSize(change_rast, minblocks = min_blocks)
      
      for (i in 1:bs$n) {
        cat("  Processing chunk", i, "of", bs$n, "\n")
        # Read chunks from both rasters
        v_current = readValues(current, bs$row[i], bs$nrows[i])
        v_previous = readValues(previous, bs$row[i], bs$nrows[i])
        
        # Calculate change
        v_change = v_current - v_previous
        
        # Apply mask: keep only where vegetation existed in at least one period
        # AND where we have valid data in both periods
        # AND mask out 0 values (no change) as requested
        valid_mask = !is.na(v_current) & !is.na(v_previous)
        vegetation_mask = (v_current == 1) | (v_previous == 1)
        change_mask = v_change != 0  # Mask out 0 values (no change)
        final_mask = valid_mask & vegetation_mask & change_mask
        
        v_change[!final_mask] = NA
        
        # Write chunk back to raster
        change_rast = writeValues(change_rast, v_change, bs$row[i])
      }
      
      change_rast = writeStop(change_rast)
      
    } else {
      # Standard processing for smaller rasters or when we have enough RAM
      change = current - previous
      
      # Create mask: keep only pixels where:
      # 1. There was vegetation in at least one period
      # 2. We have valid data in both periods  
      # 3. The change is NOT zero (mask out no-change areas)
      valid_mask = !is.na(current) & !is.na(previous)
      vegetation_mask = (current == 1) | (previous == 1)
      change_mask = change != 0  # This masks out 0 values as requested
      final_mask = valid_mask & vegetation_mask & change_mask
      
      # Apply mask - only keep areas with actual change (±1)
      change_masked = change * final_mask
      change_masked[!final_mask] = NA
      
      change_rast = change_masked
    }
    
    return(change_rast)
  }
  
  # Calculate changes for each time period
  cat("Calculating vegetation changes...\n")
  twoW_change = calculate_change_efficient(target_layer, twoW_layer, "2-week")
  oneM_change = calculate_change_efficient(target_layer, oneM_layer, "1-month")
  oneY_change = calculate_change_efficient(target_layer, oneY_layer, "1-year")
  fiveY_change = calculate_change_efficient(target_layer, fiveY_layer, "5-year")
  
  # Set names for the output rasters
  names(twoW_change) = paste0("twoW_change_", targetDATE)
  names(oneM_change) = paste0("oneM_change_", targetDATE)
  names(oneY_change) = paste0("oneY_change_", targetDATE)
  names(fiveY_change) = paste0("fiveY_change_", targetDATE)
  
  # Clear intermediate objects to free memory
  rm(target_layer, twoW_layer, oneM_layer, oneY_layer, fiveY_layer)
  gc()
  
  # Create output list
  result_list = list(
    twoW = twoW_change,
    oneM = oneM_change,
    oneY = oneY_change,
    fiveY = fiveY_change,
    metadata = list(
      target_date = target_date,
      target_name = sorted_names[target_idx],
      comparison_dates = list(
        twoW = twoW_info$date,
        oneM = oneM_info$date,
        oneY = oneY_info$date,
        fiveY = fiveY_info$date
      ),
      comparison_names = list(
        twoW = twoW_info$name,
        oneM = oneM_info$name,
        oneY = oneY_info$name,
        fiveY = fiveY_info$name
      ),
      date_differences = list(
        twoW = twoW_info$diff_days,
        oneM = oneM_info$diff_days,
        oneY = oneY_info$diff_days,
        fiveY = fiveY_info$diff_days
      ),
      processing_mode = processing_option,
      memory_fraction_used = switch(processing_option,
                                    "lowmem" = 0.3,
                                    "medium" = 0.5,
                                    "high" = 0.7,
                                    "aggressive" = 0.9,
                                    "auto" = "system_default"),
      note = "0 values (no change) are masked out as requested - only ±1 values remain"
    )
  )
  
  cat("Vegetation change calculation completed successfully!\n")
  cat("Memory mode used:", processing_option, "\n")
  cat("Note: 0 values (no change) have been masked out - only vegetation gain (+1) and loss (-1) remain\n")
  return(result_list)
}



#' Disk-based vegetation change analysis
#'
#' Memory-efficient version that processes large rasters in chunks and saves
#' results directly to disk. Ideal for very large datasets or limited RAM.
#'
#' @param inputRAST A SpatRaster object with vegetation data
#' @param targetDATE Target date for analysis obtained from names of input raster stack
#' @param temp_dir Temporary directory for processing files
#' @param auto_optimize Automatically optimize memory settings (recommended)
#' @param custom_memfrac Custom memory fraction (0-1)
#' @param custom_chunk_size Custom chunk size in pixels
#'
#' @return A list with vegetation change rasters and processing metadata
#'
#' @examples
#' \dontrun{
#' # For large datasets
#' changes <- CHUNKWISE_get_VC_TOdisk(
#'   large_veg_data, 
#'   "2023-06-15",
#'   temp_dir = "D:/temp_processing"
#' )
#' }
#'
#' @export
CHUNKWISE_get_VC_TOdisk <- function(inputRAST, 
                                    targetDATE, 
                                    temp_dir = "D:/temp_processing", 
                                    auto_optimize = TRUE,
                                    custom_memfrac = NULL,
                                    custom_chunk_size = NULL) {
  
  # Create temp directory
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Get optimal settings if auto-optimize is TRUE
  if(auto_optimize) {
    optimal_settings = CHUNKWISE_optimal_memfrac()
    mem_fraction = optimal_settings$memfrac
    chunk_size = optimal_settings$chunk_size
    cat("Auto-optimized settings applied\n")
  } else {
    # Use custom settings
    mem_fraction = ifelse(is.null(custom_memfrac), 0.6, custom_memfrac)
    chunk_size = ifelse(is.null(custom_chunk_size), 800, custom_chunk_size)
  }
  
  # Set memory options
  terraOptions(memfrac = mem_fraction, tempdir = temp_dir, tolerance = 0.1)
  
  cat("Disk-based chunk-wise processing: Using", chunk_size, "x", chunk_size, "pixel chunks\n")
  cat("Memory fraction:", mem_fraction, "(", mem_fraction*100, "% RAM)\n")
  cat("Temporary directory:", temp_dir, "\n")
  
  # Date processing
  layer_names = names(inputRAST)
  dates = as.Date(layer_names)
  
  # Remove duplicate dates by keeping the first occurrence
  unique_indices = !duplicated(dates)
  unique_rast = inputRAST[[unique_indices]]
  unique_dates = dates[unique_indices]
  unique_names = layer_names[unique_indices]
  
  # Sort by date to ensure chronological order
  date_order = order(unique_dates)
  sorted_rast = unique_rast[[date_order]]
  sorted_dates = unique_dates[date_order]
  sorted_names = unique_names[date_order]
  
  # Find the target date
  target_date = as.Date(targetDATE)
  target_idx = which(sorted_dates == target_date)
  
  if(length(target_idx) == 0) {
    stop("Target date ", targetDATE, " not found in raster layers")
  }
  
  cat("Target layer:", sorted_names[target_idx], "\n")
  
  # Function to find the closest available date with seasonal consideration
  find_seasonal_comparison = function(base_date, years_back = 0, approximate_days = NULL) {
    if (!is.null(approximate_days)) {
      # For short-term comparisons (weeks, months)
      target_lookup_date = base_date - approximate_days
    } else {
      # For year-based comparisons - try to match the same seasonal period
      target_lookup_date = base_date - years(years_back)
    }
    
    # Find the closest date to the lookup date within a reasonable window
    date_diffs = abs(sorted_dates - target_lookup_date)
    
    # For year comparisons, prioritize dates within ±30 days of the target seasonal period
    if (years_back > 0) {
      seasonal_window = 30  # days
      in_season = abs(as.numeric(sorted_dates - target_lookup_date)) <= seasonal_window
      if (any(in_season)) {
        date_diffs[!in_season] = Inf  # Penalize dates outside seasonal window
      }
    }
    
    closest_idx = which.min(date_diffs)
    
    return(list(
      idx = closest_idx, 
      date = sorted_dates[closest_idx],
      diff_days = as.numeric(abs(sorted_dates[closest_idx] - target_lookup_date)),
      name = sorted_names[closest_idx]
    ))
  }
  
  # Find comparison layers with seasonal consideration
  cat("Finding comparison dates...\n")
  twoW_info = find_seasonal_comparison(target_date, approximate_days = 14)
  oneM_info = find_seasonal_comparison(target_date, approximate_days = 30)
  oneY_info = find_seasonal_comparison(target_date, years_back = 1)
  fiveY_info = find_seasonal_comparison(target_date, years_back = 5)
  
  cat("Comparison layers found:\n")
  cat("2 weeks before:", twoW_info$name, "(difference:", twoW_info$diff_days, "days)\n")
  cat("1 month before:", oneM_info$name, "(difference:", oneM_info$diff_days, "days)\n")
  cat("1 year before:", oneY_info$name, "(difference:", oneY_info$diff_days, "days)\n")
  cat("5 years before:", fiveY_info$name, "(difference:", fiveY_info$diff_days, "days)\n")
  
  # Get the raster layers
  cat("Loading target and comparison layers...\n")
  target_layer = sorted_rast[[target_idx]]
  twoW_layer = sorted_rast[[twoW_info$idx]]
  oneM_layer = sorted_rast[[oneM_info$idx]]
  oneY_layer = sorted_rast[[oneY_info$idx]]
  fiveY_layer = sorted_rast[[fiveY_info$idx]]
  
  # Get raster dimensions
  ncols = ncol(target_layer)
  nrows = nrow(target_layer)
  
  col_chunks = ceiling(ncols / chunk_size)
  row_chunks = ceiling(nrows / chunk_size)
  total_chunks = col_chunks * row_chunks
  
  cat("Dividing raster into", total_chunks, "chunks (", col_chunks, "x", row_chunks, ")\n")
  cat("Estimated memory per chunk: ~", 
      round((chunk_size * chunk_size * 4 * 8 * 6) / 1024 / 1024, 1), "MB\n\n")
  
  # Create output files on disk
  output_files = c()
  cat("Creating output raster files on disk...\n")
  
  # Initialize output rasters with NA values using terra's efficient method
  for(period in c("twoW", "oneM", "oneY", "fiveY")) {
    output_file = file.path(temp_dir, paste0(period, "_change_", targetDATE, ".tif"))
    
    # Create output raster with same properties as target layer
    out_rast = rast(
      ext(target_layer),
      nrows = nrows,
      ncols = ncols,
      crs = crs(target_layer)
    )
    
    # Initialize with NA values efficiently
    values(out_rast) = NA
    
    # Write to disk
    writeRaster(out_rast, filename = output_file, overwrite = TRUE, NAflag = -9999)
    output_files[period] = output_file
    cat("Created:", output_file, "\n")
  }
  
  # Function to process a single chunk and write to disk 
  process_chunk_to_disk = function(row_start, row_end, col_start, col_end, chunk_id, total_chunks) {
    # Calculate chunk extent
    chunk_ext = ext(
      xFromCol(target_layer, col_start),
      xFromCol(target_layer, col_end),
      yFromRow(target_layer, row_end),
      yFromRow(target_layer, row_start)
    )
    
    # Crop all layers to this chunk
    target_chunk = crop(target_layer, chunk_ext)
    twoW_chunk = crop(twoW_layer, chunk_ext)
    oneM_chunk = crop(oneM_layer, chunk_ext)
    oneY_chunk = crop(oneY_layer, chunk_ext)
    fiveY_chunk = crop(fiveY_layer, chunk_ext)
    
    # Function to process vegetation change for a period
    process_period_chunk = function(current, previous) {
      # Calculate change (-1, 0, 1)
      change = current - previous
      
      # Calculate stable vegetation (2)
      sum_layers = current + previous
      stable = ifel(sum_layers == 2, 2, 0)
      
      # Combine both layers
      combined = change + stable
      
      # Apply masking
      valid_mask = !is.na(current) & !is.na(previous)
      vegetation_mask = (current == 1) | (previous == 1)
      final_mask = valid_mask & vegetation_mask
      
      change_final = combined * final_mask
      change_final[!final_mask] = NA
      change_final = ifel(change_final == 0, NA, change_final)
      change_final = ifel(change_final == 2, 0, change_final)  # Stable → 0
      change_final = ifel(change_final == 3, 1, change_final)  # Edge case
      
      return(change_final)
    }
    
    # Process all periods for this chunk
    twoW_result = process_period_chunk(target_chunk, twoW_chunk)
    oneM_result = process_period_chunk(target_chunk, oneM_chunk)
    oneY_result = process_period_chunk(target_chunk, oneY_chunk)
    fiveY_result = process_period_chunk(target_chunk, fiveY_chunk)
    
    # Write each result to its respective disk file using terra's mosaic approach
    periods = list(
      twoW = list(result = twoW_result, file = output_files["twoW"]),
      oneM = list(result = oneM_result, file = output_files["oneM"]),
      oneY = list(result = oneY_result, file = output_files["oneY"]),
      fiveY = list(result = fiveY_result, file = output_files["fiveY"])
    )
    
    for(period_name in names(periods)) {
      period_data = periods[[period_name]]
      
      # Read existing output raster
      output_rast = rast(period_data$file)
      
      # Get the chunk result with proper extent and resolution
      chunk_result = period_data$result
      
      # Use terra's mosaic function to combine chunks properly
      # First, create a template with NA values for the entire raster
      template = output_rast
      values(template) = NA
      
      # Crop the template to our chunk extent
      template_chunk = crop(template, chunk_ext)
      
      # Replace the template chunk values with our computed values
      values(template_chunk) = values(chunk_result)
      
      # Mosaic the chunk into the output raster
      merged = mosaic(output_rast, template_chunk, fun = "first")
      
      # Write the merged result back
      writeRaster(merged, filename = period_data$file, overwrite = TRUE)
      
      # Clean up
      rm(template, template_chunk, merged)
    }
    
    # Clean up
    rm(target_chunk, twoW_chunk, oneM_chunk, oneY_chunk, fiveY_chunk,
       twoW_result, oneM_result, oneY_result, fiveY_result)
    gc()
    
    # Progress update
    if(chunk_id %% 10 == 0 || chunk_id == total_chunks) {
      cat("Completed chunk", chunk_id, "of", total_chunks, 
          "(", round(chunk_id/total_chunks*100, 1), "%)\n")
    }
  } 
  
  # Process all chunks
  cat("Starting chunk-wise processing...\n")
  chunk_counter = 1
  
  for (row_start in seq(1, nrows, by = chunk_size)) {
    row_end = min(row_start + chunk_size - 1, nrows)
    
    for (col_start in seq(1, ncols, by = chunk_size)) {
      col_end = min(col_start + chunk_size - 1, ncols)
      
      process_chunk_to_disk(row_start, row_end, col_start, col_end, 
                            chunk_counter, total_chunks)
      
      chunk_counter = chunk_counter + 1
    }
  }
  
  # Load the final results from disk
  cat("Loading final results from disk...\n")
  twoW_output = rast(output_files["twoW"])
  oneM_output = rast(output_files["oneM"])
  oneY_output = rast(output_files["oneY"])
  fiveY_output = rast(output_files["fiveY"])
  
  # Set proper names
  names(twoW_output) = paste0("twoW_change_", targetDATE)
  names(oneM_output) = paste0("oneM_change_", targetDATE)
  names(oneY_output) = paste0("oneY_change_", targetDATE)
  names(fiveY_output) = paste0("fiveY_change_", targetDATE)
  
  # Clean up original layers
  rm(target_layer, twoW_layer, oneM_layer, oneY_layer, fiveY_layer)
  gc()
  
  # Create result list
  result_list = list(
    twoW = twoW_output,
    oneM = oneM_output,
    oneY = oneY_output,
    fiveY = fiveY_output,
    metadata = list(
      target_date = target_date,
      target_name = sorted_names[target_idx],
      comparison_dates = list(
        twoW = twoW_info$date,
        oneM = oneM_info$date,
        oneY = oneY_info$date,
        fiveY = fiveY_info$date
      ),
      processing_info = paste("Disk-based chunk-wise processing with", 
                              chunk_size, "x", chunk_size, "pixel chunks"),
      total_chunks_processed = total_chunks,
      output_files = output_files,
      temp_dir = temp_dir
    )
  )
  
  cat("\n=== PROCESSING COMPLETED SUCCESSFULLY ===\n")
  cat("Total chunks processed:", total_chunks, "\n")
  cat("Output files saved in:", temp_dir, "\n")
  cat("Final values: -1 (loss), 0 (stable), 1 (gain), NA (non-vegetated)\n")
  
  return(result_list)
}

#' Optimize memory settings for chunk-wise processing
#'
#' Automatically determines optimal memory fraction and chunk size based on
#' available system RAM.
#'
#' @param total_ram_gb Total system RAM in GB. If NULL, automatically detected.
#'
#' @return A list containing:
#' \itemize{
#'   \item{memfrac: Recommended memory fraction}
#'   \item{chunk_size: Recommended chunk size in pixels}
#'   \item{total_ram_gb: Detected total RAM}
#' }
#'
#' @examples
#' \dontrun{
#' settings <- CHUNKWISE_optimal_memfrac()
#' cat("Use chunk size:", settings$chunk_size, "x", settings$chunk_size, "\n")
#' }
#'
#' @export
CHUNKWISE_optimal_memfrac <- function(total_ram_gb = NULL) {
  # If RAM not provided, detect it automatically
  if(is.null(total_ram_gb)) {
    if(.Platform$OS.type == "windows") {
      mem_info = system("wmic computersystem get TotalPhysicalMemory", intern = TRUE)
      total_ram_gb = as.numeric(mem_info[2]) / 1024^3
    } else {
      # For Mac/Linux
      mem_info = system("sysctl hw.memsize", intern = TRUE)
      total_ram_gb = as.numeric(gsub("hw.memsize: ", "", mem_info)) / 1024^3
    }
  }
  
  # Calculate optimal memory fraction based on total RAM
  if(total_ram_gb <= 4) {
    memfrac = 0.4   # 4GB: 40% = 1.6GB
    chunk_size = 300
  } else if(total_ram_gb <= 8) {
    memfrac = 0.5   # 8GB: 50% = 4GB
    chunk_size = 500
  } else if(total_ram_gb <= 16) {
    memfrac = 0.6   # 16GB: 60% = 9.6GB
    chunk_size = 800
  } else if(total_ram_gb <= 32) {
    memfrac = 0.7   # 32GB: 70% = 22.4GB
    chunk_size = 1200
  } else {
    memfrac = 0.8   # 64GB+: 80% 
    chunk_size = 1500
  }
  
  cat("=== OPTIMAL MEMORY SETTINGS ===\n")
  cat("Total RAM:", round(total_ram_gb, 1), "GB\n")
  cat("Recommended memory fraction:", memfrac, "(", memfrac*100, "% )\n")
  cat("Recommended chunk size:", chunk_size, "x", chunk_size, "pixels\n")
  cat("Estimated RAM usage:", round(total_ram_gb * memfrac, 1), "GB\n")
  
  return(list(memfrac = memfrac, chunk_size = chunk_size, total_ram_gb = total_ram_gb))
}

#' Monitor system memory usage during processing
#'
#' Provides detailed information about system memory usage and R session
#' memory consumption. Particularly useful for monitoring memory during
#' large raster processing operations.
#'
#' @details
#' On Windows systems, this function uses WMIC to get total system RAM,
#' free RAM, and used RAM. On all platforms, it also reports R session
#' memory usage via gc().
#'
#' @return Invisible list with memory usage statistics
#'
#' @examples
#' \dontrun{
#' # Monitor memory before and after processing
#' CHUNKWISE_memo_monitor()
#' 
#' # Process large dataset
#' changes <- CHUNKWISE_get_VC_TOdisk(large_data, "2023-06-15")
#' 
#' # Monitor memory again
#' CHUNKWISE_memo_monitor()
#' }
#'
#' @export
CHUNKWISE_memo_monitor <- function() {
  cat("=== MEMORY MONITOR ===\n")
  
  if(.Platform$OS.type == "windows") {
    # Windows memory info
    mem_info = system("wmic OS get TotalVisibleMemorySize,FreePhysicalMemory /Value", intern = TRUE)
    total_mem = as.numeric(gsub("TotalVisibleMemorySize=", "", 
                                mem_info[grep("TotalVisibleMemorySize", mem_info)])) / 1024 / 1024
    free_mem = as.numeric(gsub("FreePhysicalMemory=", "", 
                               mem_info[grep("FreePhysicalMemory", mem_info)])) / 1024 / 1024
    
    cat("Total RAM:", round(total_mem, 1), "GB\n")
    cat("Free RAM:", round(free_mem, 1), "GB\n")
    cat("Used RAM:", round(total_mem - free_mem, 1), "GB\n")
  }
  
  # R memory info using gc()
  mem_info = gc()
  total_mem = sum(mem_info[, "used"]) / 1024  # Convert to MB
  cat("R memory usage:", round(total_mem, 1), "MB\n")
}

#' Save vegetation change results
#'
#' Exports vegetation change rasters and metadata to specified directory.
#'
#' @param vc_changes Vegetation changes object from analysis functions
#' @param output_dir Output directory for saving results
#'
#' @export
save_changes <- function(vc_changes, output_dir = "D:/HealthMetric") {
  
  oudir = paste0(output_dir, "/", gsub("-", "_", Sys.Date()))
  dir.create(oudir, showWarnings = FALSE)
  
  # Save each change raster
  writeRaster(vc_changes$twoW, file.path(oudir, "twoW_VCC.tif"))
  writeRaster(vc_changes$oneM, file.path(oudir, "oneM_VCC.tif")) 
  writeRaster(vc_changes$oneY, file.path(oudir, "oneY_VCC.tif"))
  writeRaster(vc_changes$fiveY, file.path(oudir, "fiveY_VCC.tif"))
  
  # Save metadata
  saveRDS(vc_changes$metadata, file.path(oudir, "metadata_VCC.rds"))
  
  cat("Change rasters saved to:", oudir, "\n")
}

#' Extract vegetation changes to polygon features
#'
#' Calculates net vegetation change statistics for polygon features using
#' precise extraction with exactextractr.
#'
#' @param vc_changes Vegetation changes object from get_VC_optimized0()
#' @param polygons SF object with polygon features
#' @param max_polygons Maximum number of polygons to process (for testing)
#'
#' @return The input polygons with added columns for mean vegetation change
#'         at each time scale
#'
#' @examples
#' \dontrun{
#' polygons_with_changes <- extract_changes_exact(changes, study_area_polygons)
#' }
#'
#' @export
#' @importFrom exactextractr exact_extract
extract_changes_exact <- function(vc_changes, polygons, max_polygons = NULL) {
  # Option to process subset for testing
  if(!is.null(max_polygons) && nrow(polygons) > max_polygons) {
    polygons = polygons[1:max_polygons, ]
    cat("Processing subset of", max_polygons, "polygons\n")
  }
  
  # Extract target date from metadata and format it
  target_date = vc_changes$metadata$target_date
  date_suffix = gsub("-", "", as.character(target_date))
  
  # Create dynamic column names
  change_cols = c(
    paste0(date_suffix, "_twoW_change"),
    paste0(date_suffix, "_oneM_change"), 
    paste0(date_suffix, "_oneY_change"),
    paste0(date_suffix, "_fiveY_change")
  )
  
  cat("Using column prefix:", date_suffix, "\n")
  
  # Transform polygons to match raster CRS
  polygons_proj = st_transform(polygons, crs(vc_changes$twoW))
  
  # Use exactextractr for precise extraction with dynamic column names
  cat("Extracting 2-week changes...\n")
  polygons_proj[[change_cols[1]]] = exactextractr::exact_extract(vc_changes$twoW, polygons_proj, 'mean')
  
  cat("Extracting 1-month changes...\n")
  polygons_proj[[change_cols[2]]] = exactextractr::exact_extract(vc_changes$oneM, polygons_proj, 'mean')
  
  cat("Extracting 1-year changes...\n")
  polygons_proj[[change_cols[3]]] = exactextractr::exact_extract(vc_changes$oneY, polygons_proj, 'mean')
  
  cat("Extracting 5-year changes...\n")
  polygons_proj[[change_cols[4]]] = exactextractr::exact_extract(vc_changes$fiveY, polygons_proj, 'mean')
  
  # Round change columns to 3 decimal places
  polygons_proj = polygons_proj %>%
    mutate(across(all_of(change_cols), ~ ifelse(is.nan(.) | is.na(.), ., round(., 3))))
  
  # Add metadata as attribute for reference
  attr(polygons_proj, "change_columns") <- change_cols
  attr(polygons_proj, "target_date") <- target_date
  
  return(polygons_proj)
}

#' Analyze vegetation change statistics
#'
#' Generates comprehensive summary statistics and visualizations for
#' vegetation changes across polygon features.
#'
#' @param polygons_with_changes SF object from extract_changes_exact()
#' @param output_dir Optional directory to save results
#'
#' @return Invisible list with summary statistics
#'
#' @examples
#' \dontrun{
#' stats <- analyze_vegetation_changes(polygons_with_changes)
#' }
#'
#' @export 
analyze_vegetation_changes <- function(polygons_with_changes, output_dir = NULL) {
  # Auto-detect change columns (columns containing "change" in name)
  change_cols = names(polygons_with_changes)[grepl("change", names(polygons_with_changes), ignore.case = TRUE)]
  
  if(length(change_cols) == 0) {
    stop("No change columns found in the data. Make sure you used extract_changes_exact() first.")
  }
  
  cat("Detected change columns:", paste(change_cols, collapse = ", "), "\n")
  
  # Summary statistics
  changes_df = polygons_with_changes %>% 
    st_drop_geometry() %>% 
    select(all_of(change_cols))
  
  # Create results list to store statistics
  results_list = list()
  
  # Prepare output string
  output_lines = c()
  output_lines = c(output_lines, "=== VEGETATION CHANGE SUMMARY ===")
  output_lines = c(output_lines, paste("Analysis date:", Sys.Date()))
  output_lines = c(output_lines, paste("Total polygons:", nrow(polygons_with_changes)))
  output_lines = c(output_lines, paste("Polygons with change data:", sum(complete.cases(changes_df))))
  output_lines = c(output_lines, paste("Change periods analyzed:", length(change_cols)))
  output_lines = c(output_lines, "")
  
  for(col in change_cols) {
    vals = changes_df[[col]]
    vals_clean = vals[!is.na(vals)]
    
    if(length(vals_clean) > 0) {
      # Calculate statistics
      stats = list(
        period = col,
        n_polygons = length(vals_clean),
        mean = round(mean(vals_clean), 4),
        sd = round(sd(vals_clean), 4),
        min = round(min(vals_clean), 4),
        max = round(max(vals_clean), 4),
        gain_percent = round(sum(vals_clean > 0) / length(vals_clean) * 100, 2),
        loss_percent = round(sum(vals_clean < 0) / length(vals_clean) * 100, 2),
        stable_percent = round(sum(vals_clean == 0) / length(vals_clean) * 100, 2)
      )
      
      results_list[[col]] = stats
      
      # Add to output lines
      output_lines = c(output_lines, paste(col, ":"))
      output_lines = c(output_lines, paste("  Mean:", stats$mean))
      output_lines = c(output_lines, paste("  SD:", stats$sd))
      output_lines = c(output_lines, paste("  Min:", stats$min))
      output_lines = c(output_lines, paste("  Max:", stats$max))
      output_lines = c(output_lines, paste("  % Gain (> 0):", stats$gain_percent, "%"))
      output_lines = c(output_lines, paste("  % Loss (< 0):", stats$loss_percent, "%"))
      output_lines = c(output_lines, paste("  % Stable (= 0):", stats$stable_percent, "%"))
      output_lines = c(output_lines, paste("  Valid polygons:", stats$n_polygons))
      output_lines = c(output_lines, "")
    }
  }
  
  cat(paste(output_lines, collapse = "\n"), "\n")
  
  # Write to files (.txt and .csv) if output_dir is provided
  if(!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    timestamp = format(Sys.time(), "%Y%m%d_%H%M%S")
    output_file = file.path(output_dir, paste0("VCC_summary_", timestamp, ".txt"))
    
    writeLines(output_lines, output_file)
    cat("Summary statistics saved to:", output_file, "\n")
    
    csv_file = file.path(output_dir, paste0("VCC_stats_", timestamp, ".csv"))
    
    if(length(results_list) > 0) {
      stats_df = do.call(rbind, lapply(results_list, as.data.frame))
      write.csv(stats_df, csv_file, row.names = FALSE)
      cat("Detailed statistics saved to:", csv_file, "\n")
    }
  }
  
  invisible(results_list)
}

#' Plot vegetation change distributions
#'
#' Creates histogram and density plots of vegetation changes across
#' different time scales.
#'
#' @param polygons_with_changes SF object from extract_changes_exact()
#'
#' @return ggplot object with the visualization
#'
#' @examples
#' \dontrun{
#' plot_polygon_changes(polygons_with_changes)
#' }
#'
#' @export
#' @import ggplot2 patchwork dplyr tidyr
plot_polygon_changes <- function(polygons_with_changes) {
  # Auto-detect change columns
  change_cols = names(polygons_with_changes)[grepl("change", names(polygons_with_changes), ignore.case = TRUE)]
  
  if(length(change_cols) == 0) {
    stop("No change columns found in the data.")
  }
  
  # Extract target date from column names (take from first column)
  first_col = change_cols[1]
  date_part = substr(first_col, 1, 8)  # YYYYMMDD part
  target_date = paste0(
    substr(date_part, 1, 4), "-", 
    substr(date_part, 5, 6), "-", 
    substr(date_part, 7, 8)
  )
  
  # Create better period labels from column names
  changes_long = polygons_with_changes %>% 
    st_drop_geometry() %>% 
    select(all_of(change_cols)) %>%
    pivot_longer(
      cols = all_of(change_cols),
      names_to = "period",
      values_to = "change_value"
    ) %>%
    mutate(
      # Extract period type from column names
      period_type = case_when(
        grepl("twoW", period) ~ "2-Week",
        grepl("oneM", period) ~ "1-Month", 
        grepl("oneY", period) ~ "1-Year",
        grepl("fiveY", period) ~ "5-Year",
        TRUE ~ "Unknown"
      ),
      # Use clean period labels for panels
      period_label = factor(period_type,
                            levels = c("2-Week", "1-Month", "1-Year", "5-Year"))
    )
  
  # Calculate summary stats for annotations
  stats_df = changes_long %>%
    group_by(period_label) %>%
    summarise(
      mean_val = mean(change_value, na.rm = TRUE),
      median_val = median(change_value, na.rm = TRUE),
      n = sum(!is.na(change_value)),
      total_polygons = n(),
      na_count = sum(is.na(change_value)),
      .groups = 'drop'
    )
  
  pp = ggplot(changes_long, aes(x = change_value)) +
    geom_hline(yintercept = 0, linewidth = .6, color = "grey78") +
    geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "steelblue", alpha = 0.7) +
    geom_density(alpha = 0.5, fill = "orange", color = "orange", linewidth = NA) +
    geom_vline(xintercept = 0, color = "grey15", linewidth = .8) +
    geom_vline(data = stats_df, aes(xintercept = mean_val), 
               color = "firebrick3", linewidth = .8) +
    facet_wrap(~ period_label, ncol = 2) +
    labs(
      title = paste("Vegetation Change Distribution (", target_date, ")"),
      subtitle = "Black line = Zero change | Red line = Mean",
      x = "Net Change Value (-1 = Complete Loss, 1 = Complete Gain)",
      y = "Density"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      strip.text = element_text(face = "bold", size = 13),
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 13),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title.x = element_text(size = 13),
      axis.text.x = element_text(size = 12)
    )
  
  print(pp)
  
  return(pp)
}

#' Helper function to recommend memory mode based on system
#'
#' Provides guidance on selecting the appropriate memory processing option
#' based on available system RAM.
#'
#' @export
recommend_memory_mode <- function(available_ram_gb = NULL) {
  if(is.null(available_ram_gb)) {
    # Simple recommendation based on common scenarios
    cat("
MEMORY MODE GUIDE:
- 'lowmem': 30% RAM - Use when running many other applications
- 'medium':  50% RAM - Good balance for normal multi-tasking  
- 'high': 70% RAM - Better performance, some room for other apps
- 'aggressive': 90% RAM - Maximum speed, minimal other tasks
- 'auto':    Let terra decide (usually 60-70%)
- If you have 16GB RAM or less: use 'lowmem' or 'medium'
- If you have 32GB RAM: use 'medium' or 'high'  
- If you have 64GB+ RAM: use 'high' or 'aggressive'
- If unsure: start with 'medium' and monitor system performance
- If you need to run other heavy applications: use 'lowmem'
- For dedicated analysis sessions: use 'aggressive'
")
  }
}

#' Monitor system memory usage
#'
#' Displays current memory usage information (Windows only).
#'
#' @export
monitor_memory <- function() {
  if(.Platform$OS.type == "windows") {
    mem = system("wmic OS get FreePhysicalMemory /Value", intern = TRUE)
    mem = as.numeric(gsub("FreePhysicalMemory=", "", mem[grep("FreePhysicalMemory", mem)])) / 1024 / 1024
    cat("Available RAM:", round(mem, 1), "GB\n")
  } else {
    cat("Memory monitoring available on Windows only\n")
  }
}

#' Clear terra cache and free memory
#'
#' Removes temporary terra files and performs garbage collection to free memory.
#'
#' @export
clear_terra_cache <- function() {
  tmp_files = list.files(tempdir(), pattern = "spat_", full.names = TRUE)
  if(length(tmp_files) > 0) {
    file.remove(tmp_files)
    cat("Removed", length(tmp_files), "temporary terra files\n")
  }
  gc()
}


