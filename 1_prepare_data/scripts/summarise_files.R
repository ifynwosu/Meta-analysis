
metadata_dir <- "/inwosu/curated_data/Data/analysis_ready_renamed_metadata"
meta_dir_paths <- list.files(metadata_dir, full.names = T)
#meta_dir_paths

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

  # big_data <- rbind (big_data, sum_ER)
  # new_data <- qpcR:::cbind.na(sum_ER, sum_PR, sum_HER2)
  # big_data <- rbind (big_data, new_data)
}

big_race_data <- race_data |>
  pivot_wider(names_from = race, values_from = n) |>
  mutate(proportion = (Black/White))

big_ER_data <- er_data |>
  pivot_wider(names_from = ER_Status, values_from = n) |>
  mutate(proportion = (negative/positive)) |>
  filter(between(proportion, 0.4, 1.7))

big_PR_data <- pr_data %>% 
  pivot_wider(names_from = PR_Status, values_from = n) |>
  mutate(proportion = (negative/positive)) |>
  filter(between(proportion, 0.1, 5.5))

big_HER2_data <- her2_data %>% 
  pivot_wider(names_from = HER2_Status, values_from = n) |>
  mutate(proportion = (negative/positive)) |>
  filter(between(proportion, 0.4, 1.7))

# write_tsv(big_race_data, file.path("/inwosu/Meta_Analysis/Data/", "race_data.tsv"))
# write_tsv(big_ER_data, file.path("/inwosu/Meta_Analysis/Data/", "ER_data.tsv"))
# write_tsv(big_PR_data, file.path("/inwosu/Meta_Analysis/Data/", "PR_data.tsv"))
# write_tsv(big_HER2_data, file.path("/inwosu/Meta_Analysis/Data/", "HER2_data.tsv"))


# code to copy select files to folder for analysis
# use this for metadata
big_PR_data$Dataset_ID <- paste0(big_PR_data$Dataset_ID, ".tsv")

# use this for expression data
big_PR_data$Dataset_ID <- paste0(big_PR_data$Dataset_ID, ".tsv.gz")

# save data set id to vector
names_file <- big_PR_data$Dataset_ID |> tibble()

old_dir <- "/inwosu/curated_data/Data/analysis_ready_renamed_metadata"
new_dir <- "/inwosu/Meta_Analysis/Data/PR_metadata"

setwd(old_dir)
for(i in names_file) {
  file.copy(i, new_dir)
  }
