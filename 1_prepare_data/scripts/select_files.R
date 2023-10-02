
# create relevant directories
race_metadata_path <- "/Data/race_metadata/"
if (!dir.exists(race_metadata_path)) {
  dir.create(race_metadata_path, recursive = TRUE)
}

race_expression_data_path <- "/Data/race_expression_data/"
if (!dir.exists(race_expression_data_path)) {
  dir.create(race_expression_data_path, recursive = TRUE)
}

ER_metadata_path <- "/Data/ER_metadata/"
if (!dir.exists(ER_metadata_path)) {
  dir.create(ER_metadata_path, recursive = TRUE)
}

ER_expression_data_path <- "/Data/ER_expression_data/"
if (!dir.exists(ER_expression_data_path)) {
  dir.create(ER_expression_data_path, recursive = TRUE)
}

PR_metadata_path <- "/Data/PR_metadata/"
if (!dir.exists(PR_metadata_path)) {
  dir.create(PR_metadata_path, recursive = TRUE)
}

PR_expression_data_path <- "/Data/PR_expression_data/"
if (!dir.exists(PR_expression_data_path)) {
  dir.create(PR_expression_data_path, recursive = TRUE)
}

HER2_metadata_path <- "/Data/HER2_metadata/"
if (!dir.exists(HER2_metadata_path)) {
  dir.create(HER2_metadata_path, recursive = TRUE)
}

HER2_expression_data_path <- "/Data/HER2_expression_data/"
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
  sum_ER <- data.frame(Dataset_ID = character(0), ER_Status = character(0), n = numeric(0))
  sum_PR <- data.frame(Dataset_ID = character(0), PR_Status = character(0), n = numeric(0))
  sum_HER2 <- data.frame(Dataset_ID = character(0), HER2_Status = character(0), n = numeric(0))
  
  file = meta_dir_paths[i]
  gseID <- file %>% basename() %>% file_path_sans_ext()
  
  metadata <- read_tsv(file) %>%
    mutate(across(everything(), as.character))
  
  if ("race" %in% names(metadata)) {
    sum_race <- metadata %>% group_by(race) %>% tally() %>% mutate(Dataset_ID = gseID) %>% relocate(Dataset_ID)
  }
  race_data <- rbind(race_data, sum_race)
  
  if ("ER_Status" %in% names(metadata)) {
    sum_ER <- metadata %>% group_by(ER_Status) %>% tally() %>% mutate(Dataset_ID = gseID) %>% relocate(Dataset_ID)
  }
  er_data <- rbind(er_data, sum_ER)
  
  if ("PR_Status" %in% names(metadata)) {
    sum_PR <- metadata %>% group_by(PR_Status) %>% tally() %>% mutate(Dataset_ID = gseID) %>% relocate(Dataset_ID)
  }
  pr_data <- rbind(pr_data, sum_PR)
  
  if ("HER2_Status" %in% names(metadata)) {
    sum_HER2 <- metadata %>% group_by(HER2_Status) %>% tally() %>% mutate(Dataset_ID = gseID) %>% relocate(Dataset_ID)
  }
  her2_data <- rbind(her2_data, sum_HER2)
}

# function doesn't work perfectly yet
# filter_data <- function(dataset, col_names, col_values, paramA, paramB) {
#   dataset |>
#   pivot_wider(names_from = col_names, values_from = col_values) |>
#     mutate(proportion = (paramA/paramB)) |>
#     arrange(proportion) |>
#     mutate(across(proportion, round, 5)) |>
#     filter(between(proportion, 0.05, 20)) |>
#     filter(paramA >= 10, paramB >= 10)
# }

# big_race_data <- filter_data(race_data, race, n, Black, White)
# big_ER_data <- filter_data(er_data, race, n, negative, positive)
# big_PR_data <- filter_data(pr_data, race, n, negative, positive)
# big_HER2_data <- filter_data(her2_data, race, n, negative, positive)

big_race_data <- race_data |>
  pivot_wider(names_from = race, values_from = n) |>
  mutate(proportion = (Black/White)) |>
  arrange(proportion) |>
  mutate(across(proportion, round, 5)) |>
  filter(between(proportion, 0.05, 20)) |>
  filter(Black >= 10, White >= 10) |>
  filter(!Dataset_ID == "GSE62944_Normal")

big_ER_data <- er_data |>
  pivot_wider(names_from = ER_Status, values_from = n) |>
  mutate(proportion = (negative/positive)) |>
  arrange(proportion) |>
  mutate(across(proportion, round, 5)) |>
  filter(between(proportion, 0.05, 20)) |>
  filter(negative >= 10, positive >= 10) |>
  filter(!Dataset_ID == "GSE62944_Normal")

big_PR_data <- pr_data %>% 
  pivot_wider(names_from = PR_Status, values_from = n) |>
  mutate(proportion = (negative/positive)) |>
  arrange(proportion) |>
  mutate(across(proportion, round, 5)) |>
  filter(between(proportion, 0.05, 20)) |>
  filter(negative >= 10, positive >= 10) |>
  filter(!Dataset_ID == "GSE62944_Normal")

big_HER2_data <- her2_data %>% 
  pivot_wider(names_from = HER2_Status, values_from = n) |>
  mutate(proportion = (negative/positive)) |>
  arrange(proportion) |>
  mutate(across(proportion, round, 5)) |>
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
big_race_data$Dataset_ID <- paste0(big_race_data$Dataset_ID, ".tsv") # add ".tsv" to dataset id
names_file <- big_race_data$Dataset_ID |> tibble()
copy_data(metadata_dir, race_metadata_path)

# expression data
big_race_data$Dataset_ID <- paste0(big_race_data$Dataset_ID, ".gz") # add ".tsv" to dataset id
names_file <- big_race_data$Dataset_ID |> tibble()
copy_data(expr_dir, race_expression_data_path)

### ER status
# metadata
big_ER_data$Dataset_ID <- paste0(big_ER_data$Dataset_ID, ".tsv") # add ".tsv" to dataset id
names_file <- big_ER_data$Dataset_ID |> tibble()
copy_data(metadata_dir, ER_metadata_path)

# expression data
big_ER_data$Dataset_ID <- paste0(big_ER_data$Dataset_ID, ".gz") # add ".tsv" to dataset id
names_file <- big_ER_data$Dataset_ID |> tibble()
copy_data(expr_dir, ER_expression_data_path)

### PR status
# metadata
big_PR_data$Dataset_ID <- paste0(big_PR_data$Dataset_ID, ".tsv") # add ".tsv" to dataset id
names_file <- big_PR_data$Dataset_ID |> tibble()
copy_data(metadata_dir, PR_metadata_path)

# expression data
big_PR_data$Dataset_ID <- paste0(big_PR_data$Dataset_ID, ".gz") # add ".tsv" to dataset id
names_file <- big_PR_data$Dataset_ID |> tibble()
copy_data(expr_dir, PR_expression_data_path)

### HER2 status
# metadata
big_HER2_data$Dataset_ID <- paste0(big_HER2_data$Dataset_ID, ".tsv") # add ".tsv" to dataset id
names_file <- big_HER2_data$Dataset_ID |> tibble()
copy_data(metadata_dir, HER2_metadata_path)

# expression data
big_HER2_data$Dataset_ID <- paste0(big_HER2_data$Dataset_ID, ".gz") # add ".tsv" to dataset id
names_file <- big_HER2_data$Dataset_ID |> tibble()
copy_data(expr_dir, HER2_expression_data_path)


#check data labels from beginning to the predicted data (GEO to predicted df)

