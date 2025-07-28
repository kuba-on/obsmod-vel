########################################

# The code below processes observational and modelled ice velocity data in a loop. The processing includes:
# 1. Data preparation: 
  # - reading glacier folders and mapping glacier names.

# 2. Observational data processing: 
  ### process_data function: 
  # - loads observational data, 
  # - filters data by date range, satellites, and time seperation between acquisition dates, 
  # - aggregates data by mean daily velocity.
  ### remove_outliers function:
  # - removes outliers from the data using IQR criterion.
  ### calculate_anomalies function:
  # - calculates anomalies using weighted means by monthly means and counts,
  # - normalises anomalies to the range [-1, 1].
  ### plot_boxplots function:
  # - plots boxplots of velocity data for each year and elevation band for initial visual inspection.
  ### create_yearly_subsets function:
  # - creates yearly subsets of the data, 
  # - normalises velocity for each yearly subset.

# 3. Model data processing:
  # - gets model data corresponding to observational data for each glacier,
  ### process_model_data function:
  # - loads model data, 
  # - filters by date range,
  # - interpolates data to daily resolution using natural cubic spline,
  # - distinguishes original modelled and interpolated data.
  ### calculate_model_anomalies function:
  # - calculates anomalies for model data and normalises them.
  ### create_yearly_subsets_model function:
  # - creates yearly subsets of the data,
  # - normalises velocity for each yearly subset.

# 4. Residuals calculation:
  ### calculate_residuals function:
  # - merges observational and modelled data by date,
  # - calculates residual anomalies between observed and modeled anomalies,
  # - normalises residual anomalies to the range [-1, 1].

# 5. Spline fitting:
  ### fit_spline_to_anomaly function:
  # - fits a smooth spline to the observed anomaly data for each elevation—yearly subset with a specified global degree of freedom,
  # - calculates residuals, standard deviation, and confidence intervals.
  ### fit_spline_to_anomaly_model function:
  # - fits a smooth spline to the modelled anomaly data for each elevation—yearly subset with a specified global degree of freedom,
  # - calculates residuals, standard deviation, and confidence intervals.
  ### fit_spline_to_residuals function:
  # - fits a smooth spline to the residual anomaly data for each elevation—yearly subset with a specified global degree of freedom,
  # - calculates residuals, standard deviation, and confidence intervals.

# 6. Calculating drivers of seasonality (ratioing):
  ### calculate_ratio function:
  # - assigns ratio between absolute areas bounded by spline functions for each elevation—yearly subset,
  # - returns ratio data centred at 1, estimating the dominant driver of seasonality in ice velocity (frontal vs fricitonal).
  
# 7. Graphing:
  # - a: plots the spline fits for observed and modelled anomalies for each elevation—yearly subset as a subplot,
  # - b: plots the spline fits for observed, modelled and residual anomalies for each elevation—yearly subset as a subplot,
  # - c: plots the ratio data for each elevation—yearly subset as a heatmap and binary heatmap (frontal or frictional).

#########################################





##### 1. DATA PREPERATION #####

# Parent directory containing all glacier folders (Gl1, Gl2, ..., GlN)
parent_directory <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2/Inputs/Observations/v2"

# Get a list of all folders in the parent directory matching the pattern "Gl<number>"
glacier_folders <- list.dirs(parent_directory, full.names = TRUE, recursive = FALSE)
glacier_folders <- glacier_folders[grepl("Gl\\d+$", glacier_folders)]

# Loop through each glacier folder
for (obs_data_path in glacier_folders) {
  # Set working directory to the current glacier folder
  setwd(obs_data_path)
  
  # Get the folder name (e.g., "Gl1")
  folder_name <- basename(obs_data_path)
  
  # Path to the CSV file that contains the mapping of ID to glacier names
  csv_path <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2/glacier_name_jess.csv"
  
  # Read the CSV file
  glacier_info <- read.csv(csv_path)
  
  # Find the glacier name corresponding to the folder (ID)
  glacier_name <- glacier_info$name[glacier_info$ID == folder_name]
  
  # Default to a generic name if no match is found
  if (length(glacier_name) == 0) {
    glacier_name <- "Unknown Glacier"
  }
  
  # Function to normalize a vector using min-max normalization
  normalize_min_max <- function(x) {
    return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  }
  
  ##### 2. OBSERVATIONAL DATA PROCESSING #####
  
  # Function to create yearly subsets
  create_yearly_subsets <- function(data, file_number) {
    years <- unique(data$year)
    yearly_data <- list()
    
    for (yr in years) {
      # Subset data for the year
      subset_data <- subset(data, year == yr)
      
      # If the subset is empty, skip
      if (nrow(subset_data) == 0) next
      
      # Normalize velocity for the subset
      subset_data$norm_vel <- normalize_min_max(subset_data$vel)
      
      # Add the year and file number to the subset name
      subset_name <- paste0("year", yr, "_", file_number)  # Using file_number for output naming
      yearly_data[[subset_name]] <- subset_data
    }
    
    return(yearly_data)
  }
  
  # Function to load, filter, and aggregate data from each file
  process_data <- function(file_name) {
    # Load the data
    data <- read.csv(file_name)
    
    # Parse date and filter the dataset by date range
    data$mid_date <- as.Date(data$mid_date, format = '%Y-%m-%d %H:%M:%S')
    
    # Set the date range (adjust this as needed for your case)
    date_start <- "2016-01-01"
    date_end <- "2021-12-31"
    
    # Filter data by date range
    data_filtered <- subset(data, mid_date >= as.Date(date_start) & mid_date <= as.Date(date_end))
    
    # Check if filtered data is empty
    if (nrow(data_filtered) == 0) {
      cat("No data for file:", file_name, "\n")
      return(NULL)
    }
    
    # Subset based on satellite and temporal conditions
    sub_data <- subset(data_filtered, 
                       (((satellite == "1A" | satellite == "1B") & 
                           (as.numeric(date_dt..days.) > 8 & as.numeric(date_dt..days.) <= 30) & 
                           ((mid_date < as.Date("2016-04-01")) | 
                              (mid_date >= as.Date("2016-10-01") & mid_date < as.Date("2017-04-01")) | 
                              (mid_date >= as.Date("2017-10-01") & mid_date < as.Date("2018-04-01")) |
                              (mid_date >= as.Date("2018-10-01") & mid_date < as.Date("2019-04-01")) |
                              (mid_date >= as.Date("2019-10-01") & mid_date < as.Date("2020-04-01")) |
                              (mid_date >= as.Date("2020-10-01") & mid_date < as.Date("2021-04-01")) |
                              (mid_date >= as.Date("2021-10-01")))) |
                          ((satellite == "2A" | satellite == "2B") & 
                             (as.numeric(date_dt..days.) > 8 & as.numeric(date_dt..days.) <= 60))), 
                       select = c("mid_date", "v..m.yr."))
    
    # Check if sub_data is empty
    if (nrow(sub_data) == 0) {
      cat("No valid subset data for file:", file_name, "\n")
      return(NULL)
    }
    
    # Convert velocity to numeric and aggregate by mean daily velocity
    sub_data$v..m.yr. <- as.numeric(sub_data$v..m.yr.)
    agg_data <- aggregate(sub_data$v..m.yr., by = list(sub_data$mid_date), FUN = mean)
    colnames(agg_data) <- c("date", "vel")
    
    # Extract year from the date and add a separate column
    agg_data$year <- as.numeric(format(agg_data$date, "%Y"))
    
    # Get the file number from the file name
    # Extracting the last number from the filename using a corrected regex
    file_number <- as.numeric(gsub("Gl\\d+_(\\d+)\\_itslive_v2_comb.csv", "\\1", file_name))  # Correctly capture the last number
    
    # Create yearly subsets
    yearly_subsets <- create_yearly_subsets(agg_data, file_number)
    
    return(yearly_subsets)
  }
  
  
  # Function to plot boxplots for each year's velocity data
  plot_boxplots <- function(yearly_data) {
    # Set up the plotting area for boxplots: 3 rows, 15 columns (45 subplots)
    par(mfrow = c(6, 15), mar = c(2, 2, 0.5, 0.2), oma = c(5, 6, 3, 0))
    
    # Define years and elevation bands to plot
    years <- c(2016, 2017, 2018, 2019, 2020, 2021)
    elevations <- 1:15  # 15 possible elevation bands
    
    # Loop through each year and elevation band to plot the data
    for (year in years) {
      for (elevation in elevations) {
        subset_name <- paste("year", year, "_", elevation, sep = "")
        
        if (!is.null(yearly_data[[subset_name]]) && "vel" %in% colnames(yearly_data[[subset_name]])) {
          # If data exists for the specific year and elevation, plot the boxplot
          boxplot(yearly_data[[subset_name]]$vel, outline = TRUE)
        } else {
          # If no data for this year/elevation, plot an empty box
          plot.new()  # Create an empty plot
          box()  # Draw box around the empty plot
        }
        
        # Add year labels on the left of the first plot in each row (first elevation)
        if (elevation == 1) {
          mtext(paste("----", year, "----", sep=""), side = 2, line = 3, cex = 1.1)
        }
        
        # Add elevation labels below the last row (for year 2021)
        if (year == 2021) {
          mtext(paste("--", (elevation * 100), "m--", sep = ""), side = 1, line = 3, cex = 1.1)
        }
      }
    }
    
    # Global labels for the whole plot
    mtext("year of study", side = 2, line = 4, outer = TRUE, cex = 1.3)
    mtext("elevation band", side = 1, line = 3, outer = TRUE, cex = 1.3)
    
    # Dynamic glacier name based on the working directory
    mtext(glacier_name, side = 3, line = 0.5, outer = TRUE, cex = 1.5)
  }
  
  # Function to remove outliers from yearly subsets
  remove_outliers <- function(yearly_data) {
    cleaned_data <- list()
    
    for (subset_name in names(yearly_data)) {
      subset_data <- yearly_data[[subset_name]]
      
      # Identify outliers based on the vel variable
      outliers <- boxplot.stats(subset_data$vel)$out
      
      # Remove outliers
      cleaned_subset <- subset(subset_data, !vel %in% outliers)
      
      # Store cleaned data
      cleaned_data[[subset_name]] <- cleaned_subset
    }
    
    return(cleaned_data)
  }
  
  # Function to calculate and normalize anomalies for each element in cleaned_data separately
  calculate_anomalies <- function(cleaned_data) {
    cleaned_data_with_anomalies <- list()
    
    # Loop through each subset in the cleaned_data
    for (subset_name in names(cleaned_data)) {
      subset_data <- cleaned_data[[subset_name]]
      
      if (!is.null(subset_data) && nrow(subset_data) > 0) {
        # Add a "month" column to the data
        subset_data$month <- format(as.Date(subset_data$date), "%Y-%m")
        
        # Calculate monthly means and counts
        monthly_summary <- aggregate(vel ~ month, data = subset_data, 
                                     FUN = function(x) c(mean = mean(x, na.rm = TRUE), count = length(x)))
        
        # Convert the results of aggregation into a data frame
        monthly_summary <- do.call(data.frame, monthly_summary)
        colnames(monthly_summary) <- c("month", "mean_vel", "count")
        
        # Merge the monthly summary back into the original dataset
        subset_data <- merge(subset_data, monthly_summary, by = "month", all.x = TRUE)
        
        # Calculate weights as the inverse of monthly observation counts
        subset_data$weight <- 1 / subset_data$count
        
        # Compute the weighted mean of velocity
        weighted_mean <- sum(subset_data$vel * subset_data$weight, na.rm = TRUE) / 
          sum(subset_data$weight, na.rm = TRUE)
        
        # Calculate anomalies as deviation from the weighted mean
        subset_data$anomaly <- subset_data$vel - weighted_mean
        
        # Scale anomalies to the range [-1, 1]
        max_abs_anomaly <- max(abs(subset_data$anomaly), na.rm = TRUE)
        if (max_abs_anomaly > 0) {
          subset_data$norm_anomaly <- subset_data$anomaly / max_abs_anomaly
        } else {
          subset_data$norm_anomaly <- 0  # If all anomalies are zero, set normalized anomalies to 0
        }
        
        # Store the subset with anomalies and normalized anomalies
        cleaned_data_with_anomalies[[subset_name]] <- subset_data
      } else {
        # Retain empty subsets if any
        cleaned_data_with_anomalies[[subset_name]] <- subset_data
      }
    }
    
    return(cleaned_data_with_anomalies)
  }
  
  # Automatically detect all files in the directory matching the pattern
  file_names <- list.files(pattern = "^Gl\\d+_(\\d+)_itslive_v2_comb\\.csv$")
  
  # Sorting the files numerically (optional but ensures correct order)
  file_names <- file_names[order(as.numeric(gsub("\\D", "", file_names)))]
  
  # Initialize an empty list to store the processed data
  all_data <- list()
  
  # Loop through each file and process the data
  for (file in file_names) {
    # Process each file and create subsets
    processed_file_data <- process_data(file)
    
    # Only append non-NULL results
    if (!is.null(processed_file_data)) {
      all_data <- c(all_data, processed_file_data)  # Append all the subsets
    }
  }
  
  # Check if any data was processed
  if (length(all_data) > 0) {
    # Combine all processed data into a single data frame if needed
    combined_data <- do.call(rbind, lapply(all_data, function(x) x[, c("date", "vel", "year", "norm_vel")]))
    
    # Print the first few rows of the combined data to check
    head(combined_data)
    
    # Plot boxplots for all yearly subsets after processing
    plot_boxplots(all_data)
    
    # Remove outliers from all data
    cleaned_data <- remove_outliers(all_data)
  } else {
    cat("No valid data was processed. Check your files or filters.\n")
  }
  
  # Assuming cleaned_data is the result of the outlier removal
  if (length(cleaned_data) > 0) {
    # Calculate anomalies after removing outliers
    data_with_anomalies_obs <- calculate_anomalies(cleaned_data)
    
    # Print the first few rows of an example subset with anomalies
    example_subset <- names(data_with_anomalies_obs)[1]
    if (!is.null(data_with_anomalies_obs[[example_subset]])) {
      print(head(data_with_anomalies_obs[[example_subset]]))
    }
  } else {
    cat("No valid data was processed. Check your files or filters.\n")
  }
  
  
  
  
  ##### MODEL DATA PROCESSING #####
  # Set the base path for model data
  model_data_base_path <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2/Inputs/Model/v1"
  
  # Get the corresponding model data path for the current folder
  model_data_path <- file.path(model_data_base_path, folder_name)
  
  # List all CSV files in the model data folder
  model_csv_files <- list.files(model_data_path, pattern = "\\.csv$", full.names = TRUE)
  print("CSV Files Found:")
  print(model_csv_files)
  
  # Function to create yearly subsets
  create_yearly_subsets_model <- function(model_data, file_number_model) {
    years <- unique(model_data$year)
    yearly_data_model <- list()
    
    for (yr in years) {
      # Subset data for the year
      subset_model_data <- subset(model_data, year == yr)
      
      # If the subset is empty, skip
      if (nrow(subset_model_data) == 0) next
      
      # Normalize velocity for the subset
      subset_model_data$norm_vel <- normalize_min_max(subset_model_data$vel)
      
      # Add the year and file number to the subset name
      subset_name_model <- paste0("year", yr, "_", file_number_model)  # Using file_number for output naming
      yearly_data_model[[subset_name_model]] <- subset_model_data
    }
    
    return(yearly_data_model)
  }
  
  # Function to load, filter, aggregate, and interpolate model data
  process_model_data <- function(file_name_model) {
    # Load the data
    model_data <- read.csv(file_name_model)
    
    # Ensure the data has the expected structure
    print(paste("Processing file:", file_name_model))
    print(head(model_data))  # Print the first few rows for inspection
    
    # Set the date range
    date_start <- "2016-01-01"
    date_end <- "2021-08-07"
    
    # Filter data by date range
    model_data_filtered <- subset(model_data, Calendar_Date >= as.Date(date_start) & Calendar_Date <= as.Date(date_end))
    
    # Check if filtered data is empty
    if (nrow(model_data_filtered) == 0) {
      cat("No data for file:", file_name_model, "\n")
      return(NULL)
    }
    
    # Reduce data to two columns
    agg_model_data <- model_data_filtered[, c("Calendar_Date", "band_mean")]
    colnames(agg_model_data) <- c("date", "vel")
    
    # Check if aggregated data is empty
    if (nrow(agg_model_data) == 0) {
      cat("No valid subset data for file:", file_name_model, "\n")
      return(NULL)
    }
    
    # Convert date and extract year
    agg_model_data$date <- as.Date(agg_model_data$date)
    agg_model_data$year <- as.numeric(format(agg_model_data$date, "%Y"))
    
    # Determine dynamic interpolation range
    dynamic_start <- min(agg_model_data$date, na.rm = TRUE)
    dynamic_end <- max(agg_model_data$date, na.rm = TRUE)
    
    # Generate the dynamic daily sequence
    full_date_seq <- seq.Date(dynamic_start, dynamic_end, by = "day")
    
    # Interpolate to daily resolution using spline interpolation
    spline_fit <- spline(
      x = agg_model_data$date,
      y = agg_model_data$vel,
      xout = full_date_seq,
      method = "natural"
    )
    
    # Create a data frame with interpolated results
    interpolated_data <- data.frame(
      date = as.Date(spline_fit$x, origin = "1970-01-01"),
      vel = spline_fit$y
    )
    interpolated_data$year <- as.numeric(format(interpolated_data$date, "%Y"))
    
    # Add a column to indicate whether a point is original or interpolated
    interpolated_data$is_original <- interpolated_data$date %in% agg_model_data$date
    
    # Get the file number from the file name using the updated regex
    file_number_model <- as.numeric(gsub("Gl\\d+_(\\d+)_model_v1.csv", "\\1", basename(file_name_model)))
    
    # Create yearly subsets with interpolated data
    yearly_subsets_model <- create_yearly_subsets_model(interpolated_data, file_number_model)
    
    return(yearly_subsets_model)
  }
  
  # Function to calculate anomalies and normalize for cleaned model data
  calculate_model_anomalies <- function(yearly_subsets_model) {
    cleaned_data_with_anomalies <- list()
    
    for (subset_name in names(yearly_subsets_model)) {
      subset_data <- yearly_subsets_model[[subset_name]]
      
      if (!is.null(subset_data) && nrow(subset_data) > 0) {
        # Calculate anomalies grouped by year
        subset_data$anomaly <- ave(subset_data$vel, subset_data$year, FUN = function(x) x - mean(x, na.rm = TRUE))
        
        # Scale anomalies to the range [-1, 1] using max absolute value
        max_abs_anomaly <- max(abs(subset_data$anomaly), na.rm = TRUE)
        subset_data$norm_anomaly <- subset_data$anomaly / max_abs_anomaly
        
        # Store the subset with anomalies and normalized anomalies
        cleaned_data_with_anomalies[[subset_name]] <- subset_data
      } else {
        cleaned_data_with_anomalies[[subset_name]] <- subset_data  # Retain empty subsets if any
      }
    }
    
    return(cleaned_data_with_anomalies)
  }
  
  # Automatically detect all files in the directory matching the pattern
  file_names_model <- list.files(model_data_path, pattern = "^Gl\\d+_\\d+_model_v1\\.csv$", full.names = TRUE)
  
  # Sorting the files numerically (optional but ensures correct order)
  file_names_model <- file_names_model[order(as.numeric(gsub("\\D", "", file_names_model)))]
  
  # Initialize an empty list to store the processed data
  all_model_data <- list()
  
  # Loop through each file and process the data
  for (file in file_names_model) {
    processed_file_model_data <- process_model_data(file)
    
    if (!is.null(processed_file_model_data)) {
      all_model_data <- c(all_model_data, processed_file_model_data)  # Append all the subsets
    }
  }
  
  # Check if any model data was processed
  if (length(all_model_data) > 0) {
    # Calculate anomalies and normalized anomalies after removing outliers
    data_with_anomalies_model <- calculate_model_anomalies(all_model_data)
    
    # Print the first few rows of an example subset with anomalies
    example_model_subset <- names(data_with_anomalies_model)[1]
    if (!is.null(data_with_anomalies_model[[example_model_subset]])) {
      print(head(data_with_anomalies_model[[example_model_subset]]))
    }
  } else {
    cat("No valid model data was processed. Check your files or filters.\n")
  }
  ##### 4. RESIDUAL CALCULATION #####
  
  # Function to calculate residuals between observed and modeled anomalies
  calculate_residuals <- function(data_with_anomalies_obs, data_with_anomalies_model) {
    residuals_data <- list()  # Initialize a list to store residual data
    
    # Get all unique keys (e.g., year_elevation) that exist in both observed and modeled datasets
    common_keys <- intersect(names(data_with_anomalies_obs), names(data_with_anomalies_model))
    
    # Loop through each common key
    for (key in common_keys) {
      obs_data <- data_with_anomalies_obs[[key]]  # Observed data for the key
      model_data <- data_with_anomalies_model[[key]]  # Modeled data for the key
      
      # Merge observed and modeled data by date
      merged_data <- merge(obs_data, model_data, by = "date", suffixes = c("_obs", "_model"), all = TRUE)
      
      # Calculate residual anomaly where both observed and modeled data exist
      merged_data$residual_anomaly <- ifelse(!is.na(merged_data$anomaly_obs) & !is.na(merged_data$anomaly_model),
                                             merged_data$anomaly_obs - merged_data$anomaly_model,
                                             NA)  # Assign NA if either anomaly is missing
      
      # Calculate residuals where both observed and modeled data exist
      merged_data$residual <- ifelse(!is.na(merged_data$vel_obs) & !is.na(merged_data$vel_model),
                                     merged_data$vel_obs - merged_data$vel_model,
                                     NA)  # Assign NA if either anomaly is missing
      
      # Calculate the max absolute value for normalization, ignoring NA values
      max_abs_residual <- max(abs(merged_data$residual_anomaly), na.rm = TRUE)
      
      # Handle the case where max_abs_residual is NA (i.e., all values are NA)
      if (is.na(max_abs_residual) || max_abs_residual == 0) {
        merged_data$residual_norm_anomaly <- NA  # Assign NA if no valid residuals exist
      } else {
        # Normalize residual anomalies to the range [-1, 1]
        merged_data$residual_norm_anomaly <- merged_data$residual_anomaly / max_abs_residual
      }
      
      # Store the residual data
      residuals_data[[key]] <- merged_data
    }
    
    return(residuals_data)
  }
  
  # Calculate residuals between observed and modeled anomalies
  if (exists("data_with_anomalies_obs") && exists("data_with_anomalies_model")) {
    residuals_data <- calculate_residuals(data_with_anomalies_obs, data_with_anomalies_model)
    
    # Print an example of the residuals for a specific key
    example_key <- names(residuals_data)[1]
    if (!is.null(residuals_data[[example_key]])) {
      print(paste("Residuals for:", example_key))
      print(head(residuals_data[[example_key]]))
    }
  } else {
    cat("Observational or modeled data with anomalies is missing. Residuals cannot be calculated.\n")
  }
  
  ##### 5. SPLINE FITTING #####
  
  # Function to fit a spline to each elevation's anomaly data (observational)
  fit_spline_to_anomaly <- function(data_with_anomalies_obs) {
    spline_fits <- list()  # To store the spline fits for each elevation
    
    # Loop through each elevation in the dataset
    for (subset_name in names(data_with_anomalies_obs)) {
      dataset <- data_with_anomalies_obs[[subset_name]]  # Get the dataset for this elevation
      
      # Check if the dataset has enough data points for spline fitting
      if (!is.null(dataset) && nrow(dataset) > 3) {
        # Fit a smooth spline to anomaly values with degrees of freedom (df) set to 40
        fit <- smooth.spline(dataset$date, dataset$anomaly, df = 10)
        
        # Calculate residuals, sigma, and confidence intervals
        res <- (fit$yin - fit$y) / (1 - fit$lev)  # Calculate residuals
        sigma <- sqrt(var(res, na.rm = TRUE))  # Estimate standard deviation of the residuals
        
        # Calculate upper and lower bounds for confidence intervals (2 sigma)
        upper <- fit$y + 2.0 * sigma * sqrt(fit$lev)
        lower <- fit$y - 2.0 * sigma * sqrt(fit$lev)
        
        # Convert the fitted x values (dates) back to Date format
        fit_dates <- as.Date(fit$x, origin = "1970-01-01")
        
        # Store the fitted results, including confidence intervals and dates
        fit_results <- data.frame(
          date = fit_dates,
          fitted_anomaly = fit$y,
          upper = upper,
          lower = lower
        )
        
        # Store the fit results for this elevation
        spline_fits[[subset_name]] <- fit_results
      } else {
        warning(paste("Not enough data to fit spline for:", subset_name))
      }
    }
    
    return(spline_fits)  # Return the list of spline fits for each elevation
  }
  
  # Example: Fitting splines to the observational data
  if (length(data_with_anomalies_obs) > 0) {
    spline_fits_obs <- fit_spline_to_anomaly(data_with_anomalies_obs)
    
    # Print the spline fit for a specific elevation (e.g., year2018_4)
    print(spline_fits_obs[["year2018_4"]])
  }
  
  
  
  
  # Function to fit a spline to each elevation's anomaly data (modeled)
  fit_spline_to_anomaly_model <- function(data_with_anomalies_model) {
    spline_fits_model <- list()  # To store the spline fits for each elevation
    
    # Loop through each elevation in the dataset
    for (subset_name in names(data_with_anomalies_model)) {
      dataset <- data_with_anomalies_model[[subset_name]]  # Get the dataset for this elevation
      
      # Check if the dataset has enough data points for spline fitting
      if (!is.null(dataset) && nrow(dataset) > 3) {
        # Fit a smooth spline to anomaly values with degrees of freedom (df) set to 15
        fit <- smooth.spline(dataset$date, dataset$anomaly, df = 50)
        
        # Calculate residuals, sigma, and confidence intervals
        res <- (fit$yin - fit$y) / (1 - fit$lev)  # Calculate residuals
        sigma <- sqrt(var(res, na.rm = TRUE))  # Estimate standard deviation of the residuals
        
        # Calculate upper and lower bounds for confidence intervals (2 sigma)
        upper <- fit$y + 2.0 * sigma * sqrt(fit$lev)
        lower <- fit$y - 2.0 * sigma * sqrt(fit$lev)
        
        # Convert the fitted x values (dates) back to Date format
        fit_dates <- as.Date(fit$x, origin = "1970-01-01")
        
        # Store the fitted results, including confidence intervals and dates
        fit_results <- data.frame(
          date = fit_dates,
          fitted_anomaly = fit$y,
          upper = upper,
          lower = lower
        )
        
        # Store the fit results for this elevation
        spline_fits_model[[subset_name]] <- fit_results
      } else {
        warning(paste("Not enough data to fit spline for:", subset_name))
      }
    }
    
    return(spline_fits_model)  # Return the list of spline fits for each elevation
  }
  
  # Example: Fitting splines to the modeled data
  if (length(data_with_anomalies_model) > 0) {
    spline_fits_model <- fit_spline_to_anomaly_model(data_with_anomalies_model)
    
    # Print the spline fit for a specific elevation (e.g., year2018_4)
    print(spline_fits_model[["year2018_4"]])
  }
  
  # Function to fit a spline to residual data
  fit_spline_to_residuals <- function(residuals_data) {
    spline_fits <- list()  # To store the spline fits for each elevation
    
    # Loop through each elevation in the residuals dataset
    for (subset_name in names(residuals_data)) {
      dataset <- residuals_data[[subset_name]]  # Get the residual dataset for this elevation
      
      # Check if the dataset has enough data points and filter out rows with NA residuals
      if (!is.null(dataset) && nrow(dataset) > 3) {
        valid_data <- subset(dataset, !is.na(residual_anomaly))  # Remove rows with NA residuals
        
        # Only proceed if there are enough valid data points for spline fitting
        if (nrow(valid_data) > 3) {
          # Fit a smooth spline to residual anomaly values
          fit <- smooth.spline(valid_data$date, valid_data$residual_anomaly, df = 10)
          
          # Calculate residuals, sigma, and confidence intervals
          res <- (fit$yin - fit$y) / (1 - fit$lev)  # Calculate residuals
          sigma <- sqrt(var(res, na.rm = TRUE))  # Estimate standard deviation of the residuals
          
          # Calculate upper and lower bounds for confidence intervals (2 sigma)
          upper <- fit$y + 2.0 * sigma * sqrt(fit$lev)
          lower <- fit$y - 2.0 * sigma * sqrt(fit$lev)
          
          # Convert the fitted x values (dates) back to Date format
          fit_dates <- as.Date(fit$x, origin = "1970-01-01")
          
          # Store the fitted results, including confidence intervals and dates
          fit_results <- data.frame(
            date = fit_dates,
            fitted_residual = fit$y,
            upper = upper,
            lower = lower
          )
          
          # Store the fit results for this elevation
          spline_fits[[subset_name]] <- fit_results
        } else {
          warning(paste("Not enough valid data to fit spline for:", subset_name))
        }
      } else {
        warning(paste("Not enough data to fit spline for:", subset_name))
      }
    }
    
    return(spline_fits)  # Return the list of spline fits for each elevation
  }
  
  # Example: Fitting splines to the residual data
  if (exists("residuals_data") && length(residuals_data) > 0) {
    spline_fits_residuals <- fit_spline_to_residuals(residuals_data)
    
    # Print the spline fit for a specific elevation (e.g., year2018_4)
    example_key <- names(spline_fits_residuals)[1]
    if (!is.null(spline_fits_residuals[[example_key]])) {
      print(paste("Spline fit for residuals:", example_key))
      print(head(spline_fits_residuals[[example_key]]))
    }
  } else {
    cat("Residuals data is missing or empty. Splines cannot be calculated.\n")
  }


  ##### 6. CALCULATING DRIVERS OF SEASONALITY (RATIOING) #####
  
  # Define years and elevation bands
  years <- 2016:2021
  elevations <- 1:15  # Assuming 15 elevation bands
  
  # Initialize matrix for storing ratios
  ratio_matrix <- matrix(NA, nrow = length(years), ncol = length(elevations),
                         dimnames = list(as.character(years), as.character(elevations * 100)))
  
  # Compute absolute areas for each year and elevation band
  for (year in years) {
    for (elevation in elevations) {
      subset_name <- paste0("year", year, "_", elevation)
      
      
      # Ensure spline fits exist
      if (!(subset_name %in% names(spline_fits_obs)) || 
          !(subset_name %in% names(spline_fits_model)) || 
          !(subset_name %in% names(spline_fits_residuals))) {
        next  # Skip if missing
      }
      
      # Extract spline fits
      spline_obs <- spline_fits_obs[[subset_name]]
      spline_model <- spline_fits_model[[subset_name]]
      spline_residual <- spline_fits_residuals[[subset_name]]
      
      # Ensure data is non-empty
      if (nrow(spline_obs) == 0 || nrow(spline_model) == 0 || nrow(spline_residual) == 0) {
        next  # Skip empty data
      }
      
      # Extract common date range
      common_dates <- Reduce(intersect, list(spline_obs$date, spline_model$date, spline_residual$date))
      
      if (length(common_dates) == 0) next  # Skip if no overlap
      
      common_dates <- as.Date(common_dates)  # ensure Date class
      spline_obs$date <- as.Date(spline_obs$date)
      spline_model$date <- as.Date(spline_model$date)
      spline_residual$date <- as.Date(spline_residual$date)
      
      # Extract values at common dates
      obs_values <- spline_obs$fitted_anomaly[spline_obs$date %in% common_dates]
      model_values <- spline_model$fitted_anomaly[spline_model$date %in% common_dates]
      residual_values <- spline_residual$fitted_residual[spline_residual$date %in% common_dates]
      
      # Remove NA values
      valid_indices <- complete.cases(obs_values, model_values, residual_values)
      
      if (sum(valid_indices) == 0) next  # Skip if all values are NA
      
      # Compute absolute areas
      area_front <- sum(abs(obs_values[valid_indices] - residual_values[valid_indices]), na.rm = TRUE)
      area_basal <- sum(abs(obs_values[valid_indices] - model_values[valid_indices]), na.rm = TRUE)
      
      # Assign ratio, avoiding division by zero
      if (area_basal == 0) {
        ratio_matrix[as.character(year), as.character(elevation * 100)] <- NA
      } else {
        ratio_matrix[as.character(year), as.character(elevation * 100)] <- area_front / area_basal
      }
    }
  }
  
  # Save as CSV
  output_dir <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2/Outputs/Ratios"
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  csv_file_path <- file.path(output_dir, paste0(folder_name, "_ratio.csv"))
  write.csv(ratio_matrix, file = csv_file_path, row.names = TRUE)
  cat("Ratio matrix saved at:", csv_file_path, "\n")
  
  
  
  # Define the subset name for year 2017, elevation 1
  subset_name <- "year2017_1"
  
  # Retrieve spline fits for the given year and elevation
  spline_obs <- spline_fits_obs[[subset_name]]
  spline_model <- spline_fits_model[[subset_name]]
  spline_residual <- spline_fits_residuals[[subset_name]]
  
  # Find common date range where all three splines exist
  common_dates <- Reduce(intersect, list(spline_obs$date, spline_model$date, spline_residual$date))
  common_dates <- as.Date(common_dates, format = '%Y-%m-%d', origin = "1970-01-01")
  
  # Extract values at common dates
  obs_values <- spline_obs$fitted_anomaly[spline_obs$date %in% common_dates]
  model_values <- spline_model$fitted_anomaly[spline_model$date %in% common_dates]
  residual_values <- spline_residual$fitted_residual[spline_residual$date %in% common_dates]
  
  # Print the values directly in RStudio
  print(data.frame(
    date = common_dates,
    obs_spline = obs_values,
    model_spline = model_values,
    residual_spline = residual_values
  ))


  
  ##### 7a. GRAPH OF OBS_MODEL ##### 
  
  # Define the path to save the PDF
  pdf_file_path <- file.path("/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2/Outputs/Obs_Model",
                             paste0(folder_name, "_om.pdf"))
  
  # Open a PDF device to save the plot (16x5 inches, landscape orientation)
  pdf(pdf_file_path, width = 18, height = 10)
  
  # Set up the plotting area: 6 rows, 15 columns
  par(mfrow = c(6, 15), mar = c(0.2, 0.5, 3, 0.2), oma = c(6, 6, 3, 0))  # Added extra space in oma for top text
  
  # Loop through each year (2016 to 2021)
  for (year in 2016:2021) {
    # Loop through each elevation band (1 to 15)
    for (elevation in 1:15) {
      # Retrieve spline fits and anomaly data
      spline_fit_obs <- spline_fits_obs[[paste0("year", year, "_", elevation)]]
      spline_fit_model <- spline_fits_model[[paste0("year", year, "_", elevation)]]
      subset_data_obs <- data_with_anomalies_obs[[paste0("year", year, "_", elevation)]]
      subset_data_model <- data_with_anomalies_model[[paste0("year", year, "_", elevation)]]
      
      # Ensure anomalies are numeric
      obs_anomalies <- as.numeric(subset_data_obs$anomaly)
      model_anomalies <- as.numeric(subset_data_model$anomaly)
      
      # Remove NA values before calculating the max absolute anomaly
      all_anomalies <- c(obs_anomalies, model_anomalies)
      all_anomalies <- all_anomalies[!is.na(all_anomalies)]
      
      # Calculate the maximum absolute value of anomalies
      if (length(all_anomalies) > 0) {
        y_max_absolute <- max(abs(all_anomalies), na.rm = TRUE)
      } else {
        y_max_absolute <- 0  # Default if no valid anomalies are found
      }
      
      # Set y-axis limits symmetrically around 0
      y_min <- -y_max_absolute
      y_max <- y_max_absolute
      
      if (!is.null(subset_data_obs) && !is.null(spline_fit_obs)) {
        # Plot the anomaly points (observed anomalies over time)
        plot(subset_data_obs$date, subset_data_obs$anomaly, xlab = "", ylab = "", main = "",
             ylim = c(y_min, y_max),  # Adjusted ylim
             xlim = as.Date(c(paste0(year, "-01-01"), paste0(year, "-12-31"))),
             yaxt = "n", xaxt = "n", pch = 1, col = adjustcolor("black", alpha.f = 0.5), cex = 0.5)
        
        # Set y-axis labels only for the first plot in each row (first elevation in the year)
        if (elevation == 1) {
          axis(2, at = c(y_min, 0, y_max), labels = c("-1.0", "0.0", "1.0"), las = 1)  # Y-axis labels
        }
        
        # Add the X-axis only to the last row (year == 2021)
        if (year == 2021) {
          axis.Date(1, at = seq.Date(as.Date(paste0(year, "-03-01")), as.Date(paste0(year, "-11-30")), by = "4 months"),
                    format = "%b", cex.axis = 1, las = 1)  # X-axis labels for March, July, November
        }
        
        # Filter spline fits for the current year
        spline_obs_filtered <- spline_fit_obs[
          spline_fit_obs$date >= as.Date(paste0(year, "-01-01")) &
            spline_fit_obs$date <= as.Date(paste0(year, "-12-31")), , drop = FALSE]
        
        # Add spline fits for observational data
        if (!is.null(spline_obs_filtered) && nrow(spline_obs_filtered) > 0) {
          polygon(c(spline_obs_filtered$date, rev(spline_obs_filtered$date)),
                  c(spline_obs_filtered$upper, rev(spline_obs_filtered$lower)),
                  col = adjustcolor("black", alpha.f = 0.15), border = NA)
          lines(spline_obs_filtered$date, spline_obs_filtered$fitted_anomaly, col = "black", lwd = 1, lty = 1)
        }
        
        # Plot modelled anomalies with differentiation between original and interpolated points
        if (!is.null(spline_fit_model) && !is.null(subset_data_model)) {
          # Separate original and interpolated points
          original_points <- subset(subset_data_model, is_original == TRUE)
          interpolated_points <- subset(subset_data_model, is_original == FALSE)
          
          # Plot interpolated points in red
          if (nrow(interpolated_points) > 0) {
            points(interpolated_points$date, interpolated_points$anomaly,
                   pch = 1, col = adjustcolor("#ffb09c", alpha.f = 0.5), cex = 0.4)
          }
          
          # Filter spline fits for the current year
          spline_model_filtered <- spline_fit_model[
            spline_fit_model$date >= as.Date(paste0(year, "-01-01")) &
              spline_fit_model$date <= as.Date(paste0(year, "-12-31")), , drop = FALSE]
          
          if (!is.null(spline_model_filtered) && nrow(spline_model_filtered) > 0) {
            polygon(c(spline_model_filtered$date, rev(spline_model_filtered$date)),
                    c(spline_model_filtered$upper, rev(spline_model_filtered$lower)),
                    col = adjustcolor("red", alpha.f = 0), border = NA)
            lines(spline_model_filtered$date, spline_model_filtered$fitted_anomaly,
                  col = adjustcolor("red", alpha.f = 0), lwd = 1, lty = 1)
          }
          
          # Plot original points in red
          if (nrow(original_points) > 0) {
            points(original_points$date, original_points$anomaly,
                   pch = 1, col = adjustcolor("red", alpha.f = 1), cex = 0.6)
          }
        }
        
        # Calculate the weighted mean for observations
        if (!is.null(subset_data_obs$vel) && !is.null(subset_data_obs$weight)) {
          weighted_mean_obs <- weighted.mean(subset_data_obs$vel, w = subset_data_obs$weight, na.rm = TRUE)
        } else {
          weighted_mean_obs <- NA
        }
        
        mean_model <- if (!is.null(subset_data_model)) mean(subset_data_model$vel, na.rm = TRUE) else NA
        normalization_factor <- max(abs(c(subset_data_obs$anomaly, subset_data_model$anomaly)), na.rm = TRUE)
        
        # Annotate weighted mean and model mean
        if (elevation == 1) {
          mtext(paste0("w mean obs:"), side = 3, line = 1.7, cex = 0.6, col = "black", adj = 0)
          mtext(paste0("mean model:"), side = 3, line = 1, cex = 0.6, col = "red", adj = 0)
          mtext(paste0("norm factor:"), side = 3, line = 0.3, cex = 0.6, col = "grey40", adj = 0)
        }
        mtext(paste0(round(weighted_mean_obs, 0)), side = 3, line = 1.7, cex = 0.6, col = "black", adj = 1)
        mtext(paste0(round(mean_model, 0)), side = 3, line = 1, cex = 0.6, col = "red", adj = 1)
        mtext(paste0(round(normalization_factor, 0)), side = 3, line = 0.3, cex = 0.6, col = "grey40", adj = 1)
        
      } else {
        # If data for this year and elevation is missing, plot a blank subplot
        plot(NULL, xlim = c(0, 1), ylim = c(-1, 1), type = "n", xaxt = "n", yaxt = "n", bty = "n")  # Blank plot
        box()  # Draw the outline around the blank plot
        if (elevation == 1) {
          axis(2, at = c(-1, 0, 1), labels = c("-1", "0", "1"), las = 1)  # Add a uniform y-axis to a blank plot
        }
        if (elevation == 1) {
          mtext(paste0("w mean obs:"), side = 3, line = 1.7, cex = 0.6, col = "black", adj = 0)
          mtext(paste0("mean model:"), side = 3, line = 1, cex = 0.6, col = "red", adj = 0)
          mtext(paste0("norm factor:"), side = 3, line = 0.3, cex = 0.6, col = "grey40", adj = 0)
        }
        if (year == 2021) {
          axis(1, at = c(0.1667, 0.5, 0.836), labels = c("Mar", "Jul", "Nov"), las = 1)  # Approximate X-axis labels
        }
      }
      
      # Optionally, add year label on the left margin next to the first plot in each row
      if (elevation == 1) {
        mtext(paste("----", year, "----", sep = ""), side = 2, line = 3, cex = 1.1)
        mtext("scaled anomaly", side = 2, line = 2.1, cex = 0.7)
      }
      
      # Optionally, add elevation labels below the last row of plots
      if (year == 2021) {
        mtext(paste("--", (elevation * 100), "m--", sep = ""), side = 1, line = 3, cex = 1.1)
      }
    }
  }
  
  # Global labels
  mtext("year of study", side = 2, line = 4, outer = TRUE, cex = 1.3)
  mtext("elevation band", side = 1, line = 4.5, outer = TRUE, cex = 1.3)
  mtext(paste0(folder_name, " - ", glacier_name), side = 3, line = 0.5, outer = TRUE, cex = 1.5)
  
  
  # Close the PDF device
  dev.off()
  
  # Notify user of the output location
  cat("Graph saved as PDF at:", pdf_file_path, "\n")
  
  
  
  
  
  
  ##### 7b. GRAPH OF OBS_MODEL_RES ##### 
  
  # Define the path to save the PDF
  pdf_file_path <- file.path("/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2/Outputs/Obs_Model_Res",
                             paste0(folder_name, "_omr.pdf"))
  
  # Open a PDF device to save the plot (16x5 inches, landscape orientation)
  pdf(pdf_file_path, width = 18, height = 10)
  
  # Set up the plotting area: 6 rows, 15 columns
  par(mfrow = c(6, 15), mar = c(0.2, 0.5, 3, 0.2), oma = c(6, 6, 3, 0))  # Added extra space in oma for top text
  
  # Loop through each year (2016 to 2021)
  for (year in 2016:2021) {
    # Loop through each elevation band (1 to 15)
    for (elevation in 1:15) {
      # Retrieve spline fits and anomaly data
      spline_fit_obs <- spline_fits_obs[[paste0("year", year, "_", elevation)]]
      spline_fit_model <- spline_fits_model[[paste0("year", year, "_", elevation)]]
      subset_data_obs <- data_with_anomalies_obs[[paste0("year", year, "_", elevation)]]
      subset_data_model <- data_with_anomalies_model[[paste0("year", year, "_", elevation)]]
      subset_residual_data <- residuals_data[[paste0("year", year, "_", elevation)]]
      spline_fit_residuals <- spline_fits_residuals[[paste0("year", year, "_", elevation)]]
      
      # Ensure anomalies are numeric
      obs_anomalies <- as.numeric(subset_data_obs$anomaly)
      model_anomalies <- as.numeric(subset_data_model$anomaly)
      residual_anomalies <- as.numeric(subset_residual_data$residual_anomaly)
      
      # Remove NA values before calculating the max absolute anomaly
      all_anomalies <- c(obs_anomalies, model_anomalies, residual_anomalies)
      all_anomalies <- all_anomalies[!is.na(all_anomalies)]
      
      # Calculate the maximum absolute value of anomalies
      if (length(all_anomalies) > 0) {
        y_max_absolute <- max(abs(all_anomalies), na.rm = TRUE)
      } else {
        y_max_absolute <- 0  # Default if no valid anomalies are found
      }
      
      # Set y-axis limits symmetrically around 0
      y_min <- -y_max_absolute
      y_max <- y_max_absolute
      
      if (!is.null(subset_data_obs) && !is.null(spline_fit_obs)) {
        # Plot the anomaly points (observed anomalies over time)
        plot(subset_data_obs$date, subset_data_obs$anomaly, xlab = "", ylab = "", main = "",
             ylim = c(y_min, y_max),  # Adjusted ylim
             xlim = as.Date(c(paste0(year, "-01-01"), paste0(year, "-12-31"))),
             yaxt = "n", xaxt = "n", pch = 1, col = adjustcolor("black", alpha.f = 0.5), cex = 0.5)
        
        # Set y-axis labels only for the first plot in each row (first elevation in the year)
        if (elevation == 1) {
          axis(2, at = c(y_min, 0, y_max), labels = c("-1.0", "0.0", "1.0"), las = 1)  # Y-axis labels
        }
        
        # Add the X-axis only to the last row (year == 2021)
        if (year == 2021) {
          axis.Date(1, at = seq.Date(as.Date(paste0(year, "-03-01")), as.Date(paste0(year, "-11-30")), by = "4 months"),
                    format = "%b", cex.axis = 1, las = 1)  # X-axis labels for March, July, November
        }
        
        # Filter spline fits for the current year
        spline_obs_filtered <- spline_fit_obs[
          spline_fit_obs$date >= as.Date(paste0(year, "-01-01")) &
            spline_fit_obs$date <= as.Date(paste0(year, "-12-31")), , drop = FALSE]
        
        # Add spline fits for observational data
        if (!is.null(spline_obs_filtered) && nrow(spline_obs_filtered) > 0) {
          polygon(c(spline_obs_filtered$date, rev(spline_obs_filtered$date)),
                  c(spline_obs_filtered$upper, rev(spline_obs_filtered$lower)),
                  col = adjustcolor("black", alpha.f = 0.15), border = NA)
          lines(spline_obs_filtered$date, spline_obs_filtered$fitted_anomaly, col = "black", lwd = 1, lty = 1)
        }
        
        # Plot modelled anomalies with differentiation between original and interpolated points
        if (!is.null(spline_fit_model) && !is.null(subset_data_model)) {
          # Separate original and interpolated points
          original_points <- subset(subset_data_model, is_original == TRUE)
          interpolated_points <- subset(subset_data_model, is_original == FALSE)
          
          # Plot interpolated points in red
          if (nrow(interpolated_points) > 0) {
            points(interpolated_points$date, interpolated_points$anomaly,
                   pch = 1, col = adjustcolor("#ffb09c", alpha.f = 0.5), cex = 0.4)
          }
          
          # Filter spline fits for the current year
          spline_model_filtered <- spline_fit_model[
            spline_fit_model$date >= as.Date(paste0(year, "-01-01")) &
              spline_fit_model$date <= as.Date(paste0(year, "-12-31")), , drop = FALSE]
          
          if (!is.null(spline_model_filtered) && nrow(spline_model_filtered) > 0) {
            polygon(c(spline_model_filtered$date, rev(spline_model_filtered$date)),
                    c(spline_model_filtered$upper, rev(spline_model_filtered$lower)),
                    col = adjustcolor("red", alpha.f = 0), border = NA)
            lines(spline_model_filtered$date, spline_model_filtered$fitted_anomaly,
                  col = adjustcolor("red", alpha.f = 0), lwd = 1, lty = 1)
          }
          
          # Plot original points in red
          if (nrow(original_points) > 0) {
            points(original_points$date, original_points$anomaly,
                   pch = 1, col = adjustcolor("red", alpha.f = 1), cex = 0.6)
          }
        }
        
        # Plot residual anomalies and spline fit if data is available
        if (!is.null(subset_residual_data) && !is.null(spline_fit_residuals) > 0) {
          # Plot the residual anomaly points
          points(subset_residual_data$date, subset_residual_data$residual_anomaly, pch = 1, col = adjustcolor("blue", alpha.f = 0.5), cex = 0.5)
        }
        
        # Add residual spline fit
        if (!is.null(spline_fit_residuals)) {
          spline_residuals_filtered <- spline_fit_residuals[
            spline_fit_residuals$date >= as.Date(paste0(year, "-01-01")) &
              spline_fit_residuals$date <= as.Date(paste0(year, "-12-31")), , drop = FALSE]
          
          if (!is.null(spline_residuals_filtered) && nrow(spline_residuals_filtered) > 0) {
            polygon(c(spline_residuals_filtered$date, rev(spline_residuals_filtered$date)),
                    c(spline_residuals_filtered$upper, rev(spline_residuals_filtered$lower)),
                    col = adjustcolor("blue", alpha.f = 0.15), border = NA)
            lines(spline_residuals_filtered$date, spline_residuals_filtered$fitted_residual, col = "blue", lwd = 1, lty = 1)
            abline(h = 0, lty = 3, lwd = 1, col = "black")  # Reference line at y = 0
          }
        }
        
        
        # Calculate the weighted mean for observations
        if (!is.null(subset_data_obs$vel) && !is.null(subset_data_obs$weight)) {
          weighted_mean_obs <- weighted.mean(subset_data_obs$vel, w = subset_data_obs$weight, na.rm = TRUE)
        } else {
          weighted_mean_obs <- NA
        }
        
        # Calculate the weighted mean for residuals
        if (!is.null(subset_residual_data$residual) && !is.null(subset_residual_data$weight)) {
          weighted_mean_res <- weighted.mean(subset_residual_data$residual, w = subset_residual_data$weight, na.rm = TRUE)
        } else {
          weighted_mean_res <- NA
        }
        
        mean_model <- if (!is.null(subset_data_model)) mean(subset_data_model$vel, na.rm = TRUE) else NA
        normalization_factor <- max(abs(c(subset_data_obs$anomaly, subset_data_model$anomaly, subset_residual_data$residual_anomaly)), na.rm = TRUE)
        
        # Annotate weighted mean and model mean
        if (elevation == 1) {
          mtext(paste0("w mean obs:"), side = 3, line = 1.7, cex = 0.6, col = "black", adj = 0)
          mtext(paste0("mean model:"), side = 3, line = 1, cex = 0.6, col = "red", adj = 0)
          #mtext(paste0("w mean res:"), side = 3, line = 1, cex = 0.6, col = "blue", adj = 0)
          mtext(paste0("norm factor:"), side = 3, line = 0.3, cex = 0.6, col = "grey40", adj = 0)
        } 
        
        mtext(paste0(round(weighted_mean_obs, 0)), side = 3, line = 1.7, cex = 0.6, col = "black", adj = 1)
        mtext(paste0(round(mean_model, 0)), side = 3, line = 1, cex = 0.6, col = "red", adj = 1)
        #mtext(paste0(round(weighted_mean_res, 0)), side = 3, line = 1, cex = 0.6, col = "blue", adj = 1)
        mtext(paste0(round(normalization_factor, 0)), side = 3, line = 0.3, cex = 0.6, col = "grey40", adj = 1)
        
      } else {
        # If data for this year and elevation is missing, plot a blank subplot
        plot(NULL, xlim = c(0, 1), ylim = c(-1, 1), type = "n", xaxt = "n", yaxt = "n", bty = "n")  # Blank plot
        box()  # Draw the outline around the blank plot
        if (elevation == 1) {
          axis(2, at = c(-1, 0, 1), labels = c("-1", "0", "1"), las = 1)  # Add a uniform y-axis to a blank plot
        }
        if (elevation == 1) {
          mtext(paste0("w mean obs:"), side = 3, line = 1.7, cex = 0.6, col = "black", adj = 0)
          mtext(paste0("mean model:"), side = 3, line = 1, cex = 0.6, col = "red", adj = 0)
          #mtext(paste0("w mean res:"), side = 3, line = 1, cex = 0.6, col = "blue", adj = 0)
          mtext(paste0("norm factor:"), side = 3, line = 0.3, cex = 0.6, col = "grey40", adj = 0)
        }
        if (year == 2021) {
          axis(1, at = c(0.1667, 0.5, 0.836), labels = c("Mar", "Jul", "Nov"), las = 1)  # Approximate X-axis labels
        }
      }
      
      # Optionally, add year label on the left margin next to the first plot in each row
      if (elevation == 1) {
        mtext(paste("----", year, "----", sep = ""), side = 2, line = 3, cex = 1.1)
        mtext("scaled anomaly", side = 2, line = 2.1, cex = 0.7)
      }
      
      # Optionally, add elevation labels below the last row of plots
      if (year == 2021) {
        mtext(paste("--", (elevation * 100), "m--", sep = ""), side = 1, line = 3, cex = 1.1)
      }
    }
  }
  
  # Global labels
  mtext("year of study", side = 2, line = 4, outer = TRUE, cex = 1.3)
  mtext("elevation band", side = 1, line = 4.5, outer = TRUE, cex = 1.3)
  mtext(paste0(folder_name, " - ", glacier_name), side = 3, line = 0.5, outer = TRUE, cex = 1.5)
  
  
  # Close the PDF device
  dev.off()
  
  # Notify user of the output location
  cat("Graph saved as PDF at:", pdf_file_path, "\n")
  
  
  
  ##### 7c. RATIO HEATMAP #####
  
  # Define file path for the PDF
  pdf_file_path <- file.path("/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2/Outputs/Heatmaps",
                             paste0(folder_name, "_hm.pdf"))
  
  # Open PDF device
  pdf(pdf_file_path, width = 9, height = 4)  # Increased width for legend space
  par(mar = c(0.2, 0.2, 0.2, 2), oma = c(4.5, 4, 1.5, 0.5))  # Adjusted margin for legend
  
  # Convert to log scale for better visualization, avoiding -Inf
  log_ratios <- log10(pmax(ratio_matrix, 10^-1))  # Clip small values at 10^-3
  
  # Define colors (blue for frontal-dominated, red for basal-dominated)
  breaks <- seq(-1, 1, length.out = 1000)  # Log scale
  color_palette <- colorRampPalette(c("#364B9A", "#EAECCC", "#A50026"))(999)
  
  # Reverse year order for plotting (2016 at top, 2021 at bottom)
  rev_years <- rev(years)
  log_ratios_reordered <- log_ratios[rev(seq_len(nrow(log_ratios))), ]  # Reverse row order
  
  # Set up plotting region with extra space for legend
  layout(matrix(c(1,2), nrow = 1), widths = c(4, 0.5))  # Main plot (4), legend (0.5)
  
  # Plot heatmap with corrected order
  image(1:length(elevations), 1:length(rev_years), t(log_ratios_reordered),
        col = color_palette, breaks = breaks, axes = FALSE,
        xlab = "", ylab = "")
  
  # Draw a black border around the entire heatmap
  rect(xleft = 0.5, xright = length(elevations) + 0.5,
       ybottom = 0.5, ytop = length(rev_years) + 0.5, 
       border = "black", lwd = 2)  # Adjust thickness as needed
  
  # Add axis labels with only tick marks (no long lines)
  axis(1, at = 1:length(elevations), labels = paste0(elevations * 100, "m"), las = 2, tck = -0.04, lwd = 0, lwd.ticks = 1, cex.axis = 1)
  axis(2, at = 1:length(rev_years), labels = rev_years, las = 1, tck = -0.04, lwd = 0, lwd.ticks = 1, cex.axis = 1)
  
  #Global labels
  mtext("year of study", side = 2, line = 3, outer = TRUE, cex = 1.1)
  mtext("elevation band", side = 1, line = 3.5, outer = TRUE, cex = 1.1)
  mtext(paste0(folder_name, " - ", glacier_name), side = 3, line = 0, outer = TRUE, cex = 1.3)
  
  ### ADDING COLOR LEGEND ###
  par(mar = c(0.2, 0, 0.2, 3))  # Adjust margin for legend panel
  image(1, seq(-1, 1, length.out = 1000), matrix(seq(-1, 1, length.out = 1000), nrow = 1),
        col = color_palette, breaks = breaks, axes = FALSE, xlab = "", ylab = "")
  
  # Define log-scale labels and positions
  log_ticks <- seq(-1, 1, by = 1)  # Positions in log10 space
  log_labels <- parse(text = paste0("10^", log_ticks))  # Log-scale labels
  
  # Add log-scale ticks and labels
  axis(4, at = log_ticks, labels = log_labels, las = 1, cex.axis = 0.8, line = 0.1)
  mtext("velocity driver", side = 4, line = 2.25, cex = 0.8)
  mtext("friction                                            front", side = 4, line = -1.3, cex = 0.8)
  
  # Draw a black border around the legend
  rect(xleft = 0.6, xright = 1.4, ybottom = min(log_ticks), ytop = max(log_ticks),
       border = "black", lwd = 2)  # Adjust thickness as needed
  
  # Close PDF device
  dev.off()
  
  cat("Heatmap saved as PDF at:", pdf_file_path, "\n")
  
  
  ##### RATIO BINARY HEATMAP #####
  
  ## Define file path for the PDF
  pdf_file_path <- file.path("/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Jess ITS_LIVE v2/Outputs/Heatmaps Binary",
                             paste0(folder_name, "_hmb.pdf"))  # Updated file name
  
  # Open PDF device
  pdf(pdf_file_path, width = 9, height = 4)  # Increased width for legend space
  par(mar = c(0.2, 0.2, 0.2, 2), oma = c(4.5, 4, 1.5, 0.5))  # Adjusted margin for legend
  
  # Convert to log scale for better visualization, avoiding -Inf
  log_ratios <- log10(pmax(ratio_matrix, 10^-1))  # Clip small values at 10^-3
  
  # Define colors (blue for frontal-dominated, red for basal-dominated)
  breaks <- seq(-1, 1, length.out = 3)  # Log scale
  color_palette <- colorRampPalette(c("#364B9A", "#A50026"))(2)
  
  # Reverse year order for plotting (2016 at top, 2021 at bottom)
  rev_years <- rev(years)
  log_ratios_reordered <- log_ratios[rev(seq_len(nrow(log_ratios))), ]  # Reverse row order
  
  # Set up plotting region with extra space for legend
  layout(matrix(c(1,2), nrow = 1), widths = c(4, 0.5))  # Main plot (4), legend (0.5)
  
  # Plot heatmap with corrected order
  image(1:length(elevations), 1:length(rev_years), t(log_ratios_reordered),
        col = color_palette, breaks = breaks, axes = FALSE,
        xlab = "", ylab = "")
  
  # Draw a black border around the entire heatmap
  rect(xleft = 0.5, xright = length(elevations) + 0.5,
       ybottom = 0.5, ytop = length(rev_years) + 0.5, 
       border = "black", lwd = 2)  # Adjust thickness as needed
  
  # Add axis labels with only tick marks (no long lines)
  axis(1, at = 1:length(elevations), labels = paste0(elevations * 100, "m"), las = 2, tck = -0.04, lwd = 0, lwd.ticks = 1, cex.axis = 1)
  axis(2, at = 1:length(rev_years), labels = rev_years, las = 1, tck = -0.04, lwd = 0, lwd.ticks = 1, cex.axis = 1)
  
  #Global labels
  mtext("year of study", side = 2, line = 3, outer = TRUE, cex = 1.1)
  mtext("elevation band", side = 1, line = 3.5, outer = TRUE, cex = 1.1)
  mtext(paste0(folder_name, " - ", glacier_name), side = 3, line = 0, outer = TRUE, cex = 1.3)
  
  ### ADDING COLOR LEGEND ###
  par(mar = c(0.2, 0, 0.2, 3))  # Adjust margin for legend panel
  image(1, seq(-1, 1, length.out = 2), matrix(seq(-1, 1, length.out = 2), nrow = 1),
        col = color_palette, breaks = breaks, axes = FALSE, xlab = "", ylab = "")
  
  # Define log-scale labels and positions
  log_ticks <- seq(-1, 1, by = 1)  # Positions in log10 space
  log_labels <- parse(text = paste0("10^", log_ticks))  # Log-scale labels
  
  # Add log-scale ticks and labels
  mtext("velocity driver", side = 4, line = 2.25, cex = 0.8)
  mtext("friction                                            front", side = 4, line = -1.3, cex = 0.8, col = "white")
  
  # Draw a black border around the legend
  rect(xleft = 0.6, xright = 1.4, ybottom = -2, ytop = 2,
       border = "black", lwd = 2)  # Adjust thickness as needed
  
  # Close PDF device
  dev.off()
  
  cat("Heatmap saved as PDF at:", pdf_file_path, "\n")
  
  cat("Processed glacier folder:", folder_name, "\n")
  
}