
metadata_dir <- "/inwosu/curated_data/Data/analysis_ready_renamed_metadata"
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

filter_data <- function(dataset, col_names, col_values, paramA, paramB) {
  dataset |>
  pivot_wider(names_from = col_names, values_from = col_names) |>
    mutate(proportion = (paramA/paramB)) |>
    arrange(proportion) |>
    mutate(across(proportion, round, 5)) |>
    filter(between(proportion, 0.05, 20)) |>
    filter(paramA >= 10, paramB >= 10)
}

big_race_data <- filter_data(race_data, race, n, Black, White)

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

# write_tsv(big_race_data, file.path("/inwosu/Meta_Analysis/Data/", "race_data.tsv"))
# write_tsv(big_ER_data, file.path("/inwosu/Meta_Analysis/Data/", "ER_data.tsv"))
# write_tsv(big_PR_data, file.path("/inwosu/Meta_Analysis/Data/", "PR_data.tsv"))
# write_tsv(big_HER2_data, file.path("/inwosu/Meta_Analysis/Data/", "HER2_data.tsv"))


# code to copy select files to folder for analysis
# use this for metadata
big_race_data$Dataset_ID <- paste0(big_race_data$Dataset_ID, ".tsv")

# save dataset id to vector
names_file <- big_race_data$Dataset_ID |> tibble()

old_dir <- "/inwosu/curated_data/Data/analysis_ready_renamed_metadata/"
new_dir <- "/inwosu/Meta_Analysis/Data/race_metadata/" 

#old_dir <- "/inwosu/Meta_Analysis/Data/" 
setwd(old_dir)
for(i in names_file) {
  file.copy(i, new_dir)
}

# use this for expression data
big_race_data$Dataset_ID <- paste0(big_race_data$Dataset_ID, ".gz")

# save data set id to vector
names_file <- big_race_data$Dataset_ID |> tibble()

old_dir <- "/inwosu/curated_data/Data/analysis_ready_expression_data/"
new_dir <- "/inwosu/Meta_Analysis/Data/race_expression_data/"

setwd(old_dir)
for(i in names_file) {
  file.copy(i, new_dir)
}

#check data labels from beginning to the predicted data (GEO to predicted df)
