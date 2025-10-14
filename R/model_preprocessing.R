########################################

# The code below preprocesses modelled velocity output from ISSM. The preprocessing involves:
# - Adding row names,
# - Converting from wide to long format,
# - Converting decimal dates to calendar dates,
# - Calculating the mean velocity across all elevation bands for each date.

########################################


install.packages("lubridate")
library(lubridate)

# Base directories
input_base_dir <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Data/Sampling_Points_Interpolated_Model"
output_base_dir_1 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Inputs/Model/vel_fric1"
output_base_dir_2 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Inputs/Model/vel_fric2"
output_base_dir_3 <- "/Users/jagon/Documents/Projects/Collabs/Jessica Badgeley/Version 3/Inputs/Model/vel_fric3"

# Function to preprocess a single file
preprocess_file <- function(file_name) {
  jab <- read.csv(file_name, header = FALSE)
  
  # Add row names
  row_names <- c(paste0("Point_", 1:(nrow(jab) - 1)), "Decimal_Date")
  jab <- cbind(row_names, jab)
  rownames(jab) <- jab[, 1]
  jab <- jab[, -1]
  jab <- rbind(jab[nrow(jab), ], jab[-nrow(jab), ])
  
  # Transpose and convert to data frame
  t_jab <- t(jab)
  t_jab <- as.data.frame(t_jab)
  
  # Ensure Decimal_Date is numeric
  t_jab$Decimal_Date <- as.numeric(as.character(t_jab$Decimal_Date))
  
  # Convert all other columns to numeric
  for (col in 2:ncol(t_jab)) {
    t_jab[, col] <- as.numeric(as.character(t_jab[, col]))
  }
  
  # Convert decimal date to formatted date
  date <- date_decimal(t_jab$Decimal_Date)
  calendar_date <- format(date, "%Y-%m-%d")
  
  # Calculate the mean of the measurements (excl. the first column)
  t_jab$band_mean <- rowMeans(t_jab[, 2:ncol(t_jab)], na.rm = TRUE)
  
  # Add the formatted date as a new column
  t_jab$Calendar_Date <- calendar_date
  
  # Reorder columns
  t_jab <- t_jab[, c("Decimal_Date", "Calendar_Date", setdiff(names(t_jab), c("Decimal_Date", "Calendar_Date")))]
  
  return(t_jab)
}

# List all glacier folders
glacier_folders <- list.dirs(input_base_dir, full.names = TRUE, recursive = FALSE)

# Batch proceessing of all glacier folders
for (input_dir in glacier_folders) {
  # Extract folder name and convert to output format (spaces instead of underscores)
  folder_name <- basename(input_dir)
  output_folder_name <- gsub("_", " ", folder_name)
  
  # Define corresponding output directories
  output_dir_1 <- file.path(output_base_dir_1, output_folder_name)
  output_dir_2 <- file.path(output_base_dir_2, output_folder_name)
  output_dir_3 <- file.path(output_base_dir_3, output_folder_name)
  
  # Create output directories if needed
  dir.create(output_dir_1, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir_2, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir_3, recursive = TRUE, showWarnings = FALSE)
  
  # List matching CSV files
  file_1_names <- list.files(input_dir, pattern = "_vel_fric1.csv$", full.names = TRUE)
  file_2_names <- list.files(input_dir, pattern = "_vel_fric2.csv$", full.names = TRUE)
  file_3_names <- list.files(input_dir, pattern = "_vel_fric3.csv$", full.names = TRUE)
  
  # Process and save file_1
  for (file_1_name in file_1_names) {
    processed_data <- preprocess_file(file_1_name)
    new_filename <- paste0(tools::file_path_sans_ext(basename(file_1_name)), "_procd.csv")
    write.csv(processed_data, file = file.path(output_dir_1, new_filename), row.names = FALSE)
  }
  
  # Process and save file_2
  for (file_2_name in file_2_names) {
    processed_data <- preprocess_file(file_2_name)
    new_filename <- paste0(tools::file_path_sans_ext(basename(file_2_name)), "_procd.csv")
    write.csv(processed_data, file = file.path(output_dir_2, new_filename), row.names = FALSE)
  }
  
  # Process and save file_3
  for (file_3_name in file_3_names) {
    processed_data <- preprocess_file(file_3_name)
    new_filename <- paste0(tools::file_path_sans_ext(basename(file_3_name)), "_procd.csv")
    write.csv(processed_data, file = file.path(output_dir_3, new_filename), row.names = FALSE)
  }
}



