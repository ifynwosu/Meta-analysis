library(tidyverse)
library(tools)

# Define the file path to the metadata directory
metadata_dir <- "/Data/renamed_metadata/"

# Create the metadata folder if it doesn't exist
if (!dir.exists(metadata_dir)) {
  dir.create(metadata_dir, recursive = TRUE)
}

meta_dir <- "/Data/prelim_metadata"
#meta_dir <- "/inwosu/curated_data/Data/prelim_metadata"
meta_dir_paths <- list.files(meta_dir, full.names = T)
#meta_dir_paths

for (i in seq_along(meta_dir_paths)) {
  file = meta_dir_paths[i]
  gseID <- file %>% basename() %>% file_path_sans_ext()
  
  metadata <- read_tsv(file) %>%
    mutate(across(everything(), as.character))
  
  if (gseID == "E_TABM_158") {
    new_metadata <- metadata %>%
      rename(race = "ethnicity") %>%
      replace_race()
  }
  
  if (gseID == "GSE86374") {
    new_metadata <- metadata %>%
      rename(race = "population") %>%
      replace_race()
  }
  
  if (gseID == "GSE19697") {
    new_metadata <- metadata %>%
      replace_race()
  }
  
  if (gseID == "GSE20194") {
    new_metadata <- metadata %>%
      replace_race()
  }
  
  if (gseID == "GSE20271") {
    new_metadata <- metadata %>%
      replace_race()
  }
  
  if (gseID == "GSE24185") {
    new_metadata <- metadata %>%
      replace_race()
  }

  if (gseID == "GSE50948") {
    new_metadata <- metadata %>%
      replace_race()
  }
  
  if (gseID == "GSE5847") {
    new_metadata <- metadata %>%
      replace_race()
  }
  
  if (gseID == "GSE62944_Normal") {
    new_metadata <- metadata %>%
      replace_race()
  }
  
  if (gseID == "GSE62944_Tumor") {
    new_metadata <- metadata %>%
      replace_race()
  }
  
  if (gseID == "GSE76275") {
    new_metadata <- metadata %>%
      replace_race()
  }
  
  write_tsv(new_metadata, file.path(
    metadata_dir,
    paste0(gseID, ".tsv")))
  
  print(paste0("Saved ", gseID, ".tsv to ", metadata_dir))
    
}

