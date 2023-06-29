# Define the file path to the data directory
data_dir <- "/Data/data"

# Create the data folder if it doesn't exist
if (!dir.exists(data_dir)) {
  dir.create(data_dir)
}

# Define the file path to the plots directory
plots_dir <- "plots" 

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results" 

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}