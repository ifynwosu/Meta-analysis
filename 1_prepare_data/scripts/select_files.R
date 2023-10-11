
# create relevant directories
race_metadata_path <- "/Data/metadata/race/"
if (!dir.exists(race_metadata_path)) {
  dir.create(race_metadata_path, recursive = TRUE)
}

race_expression_data_path <- "/Data/expression_data/race/"
if (!dir.exists(race_expression_data_path)) {
  dir.create(race_expression_data_path, recursive = TRUE)
}

ER_metadata_path <- "/Data/metadata/ER_status/"
if (!dir.exists(ER_metadata_path)) {
  dir.create(ER_metadata_path, recursive = TRUE)
}

ER_expression_data_path <- "/Data/expression_data/ER_status/"
if (!dir.exists(ER_expression_data_path)) {
  dir.create(ER_expression_data_path, recursive = TRUE)
}

PR_metadata_path <- "/Data/metadata/PR_status/"
if (!dir.exists(PR_metadata_path)) {
  dir.create(PR_metadata_path, recursive = TRUE)
}

PR_expression_data_path <- "/Data/expression_data/PR_status/"
if (!dir.exists(PR_expression_data_path)) {
  dir.create(PR_expression_data_path, recursive = TRUE)
}

HER2_metadata_path <- "/Data/metadata/HER2_status/"
if (!dir.exists(HER2_metadata_path)) {
  dir.create(HER2_metadata_path, recursive = TRUE)
}

HER2_expression_data_path <- "/Data/expression_data/HER2_status/"
if (!dir.exists(HER2_expression_data_path)) {
  dir.create(HER2_expression_data_path, recursive = TRUE)
}

# directory where expression data is stored
expr_dir <- "/Data/analysis_ready_expression_data/"

# assign directory where clean metadata is stored
metadata_dir <- "/Data/analysis_ready_renamed_metadata/"
meta_dir_paths <- list.files(metadata_dir, full.names = T)

race_data <- tibble()
er_data <- tibble()
pr_data <- tibble()
her2_data <- tibble()

for (i in seq_along(meta_dir_paths)) {
  sum_race <- data.frame(Dataset_ID = character(0), race = character(0), n = numeric(0))
  sum_ER <- data.frame(Dataset_ID = character(0), ER_status = character(0), n = numeric(0))
  sum_PR <- data.frame(Dataset_ID = character(0), PR_status = character(0), n = numeric(0))
  sum_HER2 <- data.frame(Dataset_ID = character(0), HER2_status = character(0), n = numeric(0))
  
  file = meta_dir_paths[i]
  gseID <- file %>% basename() %>% file_path_sans_ext()
  
  metadata <- read_tsv(file) %>%
    mutate(across(everything(), as.character))
  
  if ("race" %in% names(metadata)) {
    sum_race <- metadata %>% group_by(race) %>% tally() %>% mutate(Dataset_ID = gseID) %>% relocate(Dataset_ID)
  }
  race_data <- rbind(race_data, sum_race)
  
  if ("ER_status" %in% names(metadata)) {
    sum_ER <- metadata %>% group_by(ER_status) %>% tally() %>% mutate(Dataset_ID = gseID) %>% relocate(Dataset_ID)
  }
  er_data <- rbind(er_data, sum_ER)
  
  if ("PR_status" %in% names(metadata)) {
    sum_PR <- metadata %>% group_by(PR_status) %>% tally() %>% mutate(Dataset_ID = gseID) %>% relocate(Dataset_ID)
  }
  pr_data <- rbind(pr_data, sum_PR)
  
  if ("HER2_status" %in% names(metadata)) {
    sum_HER2 <- metadata %>% group_by(HER2_status) %>% tally() %>% mutate(Dataset_ID = gseID) %>% relocate(Dataset_ID)
  }
  her2_data <- rbind(her2_data, sum_HER2)
}

# save race and hormone receptor information
# big_race_data <- filter_data(race_data, race, n, Black, White)
# big_ER_data <- filter_data(er_data, race, n, negative, positive)
# big_PR_data <- filter_data(pr_data, race, n, negative, positive)
# big_HER2_data <- filter_data(her2_data, race, n, negative, positive)
# write_tsv(big_race_data, file.path("/Data/", "race_data.tsv")) (use as template)

big_race_data <- race_data |>
  pivot_wider(names_from = race, values_from = n) |>
  mutate(proportion = (Black/White)) |>
  arrange(proportion) |>
  #mutate(across(proportion, round, 5)) |>
  filter(between(proportion, 0.05, 20)) |>
  filter(Black >= 10, White >= 10) |>
  filter(!Dataset_ID == "GSE62944_Normal")

big_ER_data <- er_data |>
  pivot_wider(names_from = ER_status, values_from = n) |>
  mutate(proportion = (negative/positive)) |>
  arrange(proportion) |>
  #mutate(across(proportion, round, 5)) |>
  filter(between(proportion, 0.05, 20)) |>
  filter(negative >= 10, positive >= 10) |>
  filter(!Dataset_ID == "GSE62944_Normal")

big_PR_data <- pr_data %>% 
  pivot_wider(names_from = PR_status, values_from = n) |>
  mutate(proportion = (negative/positive)) |>
  arrange(proportion) |>
  #mutate(across(proportion, round, 5)) |>
  filter(between(proportion, 0.05, 20)) |>
  filter(negative >= 10, positive >= 10) |>
  filter(!Dataset_ID == "GSE62944_Normal")

big_HER2_data <- her2_data %>% 
  pivot_wider(names_from = HER2_status, values_from = n) |>
  mutate(proportion = (negative/positive)) |>
  arrange(proportion) |>
  #mutate(across(proportion, round, 5)) |>
  filter(between(proportion, 0.05, 20)) |>
  filter(negative >= 10, positive >= 10) |>
  filter(!Dataset_ID == "GSE62944_Normal")

# write_tsv(big_race_data, file.path("/Data/", "race_data.tsv"))
# write_tsv(big_ER_data, file.path("/Data/", "ER_data.tsv"))
# write_tsv(big_PR_data, file.path("/Data/", "PR_data.tsv"))
# write_tsv(big_HER2_data, file.path("/Data/", "HER2_data.tsv"))

# copy necessary files to appropriate folders for analysis in later steps
# function to copy relevant files
copy_data <- function(data_from, data_to) {
  for(i in names_file) {
    dataPath <- paste0(data_from, i)
    file.copy(dataPath, data_to)
  }
}

### race
# metadata
big_race_data$Dataset_ID <- paste0(big_race_data$Dataset_ID, ".tsv") # add ".tsv" extention to dataset id
names_file <- big_race_data$Dataset_ID |> tibble()
copy_data(metadata_dir, race_metadata_path)

# expression data
big_race_data$Dataset_ID <- paste0(big_race_data$Dataset_ID, ".gz") # add ".gz" extention to dataset id
names_file <- big_race_data$Dataset_ID |> tibble()
copy_data(expr_dir, race_expression_data_path)

### ER status
# metadata
big_ER_data$Dataset_ID <- paste0(big_ER_data$Dataset_ID, ".tsv") 
names_file <- big_ER_data$Dataset_ID |> tibble()
copy_data(metadata_dir, ER_metadata_path)

# expression data
big_ER_data$Dataset_ID <- paste0(big_ER_data$Dataset_ID, ".gz") 
names_file <- big_ER_data$Dataset_ID |> tibble()
copy_data(expr_dir, ER_expression_data_path)

### PR status
# metadata
big_PR_data$Dataset_ID <- paste0(big_PR_data$Dataset_ID, ".tsv") 
names_file <- big_PR_data$Dataset_ID |> tibble()
copy_data(metadata_dir, PR_metadata_path)

# expression data
big_PR_data$Dataset_ID <- paste0(big_PR_data$Dataset_ID, ".gz") 
names_file <- big_PR_data$Dataset_ID |> tibble()
copy_data(expr_dir, PR_expression_data_path)

### HER2 status
# metadata
big_HER2_data$Dataset_ID <- paste0(big_HER2_data$Dataset_ID, ".tsv")
names_file <- big_HER2_data$Dataset_ID |> tibble()
copy_data(metadata_dir, HER2_metadata_path)

# expression data
big_HER2_data$Dataset_ID <- paste0(big_HER2_data$Dataset_ID, ".gz") 
names_file <- big_HER2_data$Dataset_ID |> tibble()
copy_data(expr_dir, HER2_expression_data_path)


#check data labels from beginning to the predicted data (GEO to predicted df)

