create_directory <- function(directory_path) {
    if (!dir.exists(directory_path)) {
        dir.create(directory_path, recursive = TRUE)
    }
    return(directory_path)
}

copy_data <- function(names_file, data_from, data_to) {
    for (i in names_file) {
    # for (i in slice(names_file,-1)) {
        dataPath <- paste0(data_from, i)
        file.copy(dataPath, data_to)
    }
}

identify_trip_neg <- function(data) {
    data |>
        mutate(tri_neg_status = case_when(
            ER_status == "positive" ~ "non_tri_neg",
            PR_status == "positive" ~ "non_tri_neg",
            HER2_status == "positive" ~ "non_tri_neg",
            ER_status == "negative" & PR_status == "negative" & HER2_status == "negative" ~ "tri_neg")) |>
        select(Dataset_ID, Sample_ID, ER_status, PR_status, HER2_status, tri_neg_status, everything())
}

# function to copy relevant files to metadata and expression data to appropriate folders for analysis in later steps
process_data <- function(data, metadata_dir, metadata_path, expr_dir, expression_data_path) {    
    
     # Add ".tsv" extension to dataset id
    data$Dataset_ID <- paste0(data$Dataset_ID, ".tsv")
    names_file <- data$Dataset_ID |> as_tibble() |> filter(row_number() <= n() - 1)  
    copy_data(names_file, metadata_dir, metadata_path)
    
    # Add ".gz" extension to dataset id + .tsv from above
    data$Dataset_ID <- paste0(data$Dataset_ID, ".gz")
    names_file <- data$Dataset_ID |> as_tibble() |> filter(row_number() <= n() - 1)
    copy_data(names_file, expr_dir, expression_data_path)
}

# directory where standardized metadata is stored
metadata_dir <- "/Data/temp_metadata/"
meta_dir_paths <- list.files(metadata_dir, full.names = T)

# directory where scaled expression data is stored
expr_dir <- "/Data/temp_expression_data/"

# paths to store data for next step
metadata_path <- "/Data/temp2_metadata/"
expression_data_path <- "/Data/temp2_expr_data/"

# create different folders based on HR status
race_metadata_path <- create_directory(paste0(metadata_path, "race/"))
race_expression_data_path <- create_directory(paste0(expression_data_path, "race/")) 

ER_metadata_path <- create_directory(paste0(metadata_path, "ER_status/"))
ER_expression_data_path <- create_directory(paste0(expression_data_path, "ER_status/"))

PR_metadata_path <- create_directory(paste0(metadata_path, "PR_status/"))
PR_expression_data_path <- create_directory(paste0(expression_data_path, "PR_status/"))

HER2_metadata_path <- create_directory(paste0(metadata_path, "HER2_status/"))
HER2_expression_data_path <- create_directory(paste0(expression_data_path, "HER2_status/"))

triple_neg_metadata_path <- create_directory(paste0(metadata_path, "tri_neg_status/"))
triple_neg_expr_data_path <- create_directory(paste0(expression_data_path, "tri_neg_status/"))

race_tri_neg_metadata_path <- create_directory(paste0(metadata_path, "race_tri_neg/"))
race_tri_neg_expr_data_path <- create_directory(paste0(expression_data_path, "race_tri_neg/"))

# path to store summary data
summary_data_path <- create_directory("/Data/summary/")

#create empty tibbles to hold information
race_data <- tibble()
er_data <- tibble()
pr_data <- tibble()
her2_data <- tibble()

for (i in seq_along(meta_dir_paths)) {

    sum_race <- data.frame(Dataset_ID = character(0), race = character(0), n = numeric(0))
    sum_ER <- data.frame(Dataset_ID = character(0), ER_status = character(0), n = numeric(0))
    sum_PR <- data.frame(Dataset_ID = character(0), PR_status = character(0), n = numeric(0))
    sum_HER2 <- data.frame(Dataset_ID = character(0), HER2_status = character(0), n = numeric(0))
    
    file <- meta_dir_paths[i]
    gseID <- file |> basename() |> file_path_sans_ext()
    
    metadata <- read_tsv(file) |>
        mutate(across(everything(), as.character))
    
    if ("race" %in% names(metadata)) {
        sum_race <- metadata |> group_by(race) |> tally() |> mutate(Dataset_ID = gseID) |> relocate(Dataset_ID)
    }
    race_data <- rbind(race_data, sum_race)
    
    if ("ER_status" %in% names(metadata)) {
        sum_ER <- metadata |> group_by(ER_status) |> tally() |> mutate(Dataset_ID = gseID) |> relocate(Dataset_ID)
    }
    er_data <- rbind(er_data, sum_ER)
    
    if ("PR_status" %in% names(metadata)) {
        sum_PR <- metadata |> group_by(PR_status) |> tally() |> mutate(Dataset_ID = gseID) |> relocate(Dataset_ID)
    }
    pr_data <- rbind(pr_data, sum_PR)
    
    if ("HER2_status" %in% names(metadata)) {
        sum_HER2 <- metadata |> group_by(HER2_status) |> tally() |> mutate(Dataset_ID = gseID) |> relocate(Dataset_ID)
    }
    her2_data <- rbind(her2_data, sum_HER2)
}

big_race_data <- race_data |>
    filter(race %in% c("Black", "White")) |>
    pivot_wider(names_from = race, values_from = n) |>
    mutate(proportion = (Black/White)) |>
    arrange(proportion) |>
    filter(between(proportion, 0.05, 20)) |>
    filter(Black >= 10, White >= 10) |>
    janitor::adorn_totals()

big_ER_data <- er_data |>
    pivot_wider(names_from = ER_status, values_from = n) |>
    mutate(proportion = (negative/positive)) |>
    arrange(proportion) |>
    filter(between(proportion, 0.05, 20)) |>
    filter(negative >= 10, positive >= 10) |>
    janitor::adorn_totals()

big_PR_data <- pr_data |> 
    pivot_wider(names_from = PR_status, values_from = n) |>
    mutate(proportion = (negative/positive)) |>
    arrange(proportion) |>
    filter(between(proportion, 0.05, 20)) |>
    filter(negative >= 10, positive >= 10) |>
    janitor::adorn_totals()

big_HER2_data <- her2_data |>
    pivot_wider(names_from = HER2_status, values_from = n) |>
    mutate(proportion = (negative/positive)) |>
    arrange(proportion) |>
    #mutate(across(proportion, round, 5)) |>
    filter(between(proportion, 0.05, 20)) |>
    filter(negative >= 10, positive >= 10) |>
    #filter(!Dataset_ID == "GSE62944_Normal") |>
    janitor::adorn_totals()

## Identify datasets with all 3 HR statuses
ER <- big_ER_data$Dataset_ID |> as_tibble()
PR <- big_PR_data$Dataset_ID |> as_tibble()
HER2 <- big_HER2_data$Dataset_ID |> as_tibble()

common_HR_status <- list(ER, PR, HER2)
big_tri_neg_data <- reduce(common_HR_status, inner_join) |> 
    `colnames<-`("Dataset_ID") |> 
    filter(row_number() <= n() - 1)                                   # remove last row (contains "Totals" from janitor::adorn_totals)
    
process_data(big_race_data, metadata_dir, race_metadata_path, expr_dir, race_expression_data_path)
process_data(big_ER_data, metadata_dir, ER_metadata_path, expr_dir, ER_expression_data_path)
process_data(big_PR_data, metadata_dir, PR_metadata_path, expr_dir, PR_expression_data_path)
process_data(big_HER2_data, metadata_dir, HER2_metadata_path, expr_dir, HER2_expression_data_path)
process_data(big_tri_neg_data, metadata_dir, triple_neg_metadata_path, expr_dir, triple_neg_expr_data_path)

## identifying triple negative samples
triple_negative_files <- list.files(triple_neg_metadata_path, full.names = TRUE)
trip_neg_data <- tibble()

for (i in seq_along(triple_negative_files)) {    
    sum_trip_neg <- data.frame(Dataset_ID = character(0), tri_neg_status = character(0), n = numeric(0))    
    meta_file <- triple_negative_files[i]
    gseID <- meta_file |> 
        basename() |> 
        file_path_sans_ext()

    meta_data <- read_tsv(meta_file) |>
        mutate(across(everything(), as.character)) |>
        identify_trip_neg() |>
        filter(!is.na(tri_neg_status))

    sum_trip_neg <- meta_data |> group_by(tri_neg_status) |> tally() |> mutate(Dataset_ID = gseID) |> relocate(Dataset_ID)    
    trip_neg_data <- rbind(trip_neg_data, sum_trip_neg)

    write_tsv(meta_data, file.path(triple_neg_metadata_path,  paste0(meta_data$Dataset_ID[1], ".tsv")))
}

big_tri_neg_data <- trip_neg_data |>
    pivot_wider(names_from = tri_neg_status, values_from = n) |>
    janitor::adorn_totals()

## identify datasets with race and triple negative status
race_tri_neg_data <- inner_join(big_race_data, big_tri_neg_data)  # This is to get a list of common datasets with race and tri-neg data
race_tri_neg <- tibble()

for (i in seq_along(triple_negative_files)) {
    sum_race_tri_neg <- data.frame(Dataset_ID = character(0), tri_neg_status = character(0), n = numeric(0))
    file <- triple_negative_files[i]
    gseID <- file |> 
        basename() |> 
        file_path_sans_ext()
  
    metadata <- read_tsv(file) |>        
        mutate(across(everything(), as.character))
  
    if ("race" %in% names(metadata)) {
        sum_race_tri_neg <- metadata |> group_by(race, tri_neg_status) |> tally() |> mutate(Dataset_ID = gseID) |> relocate(Dataset_ID)
    }
    race_tri_neg <- rbind(race_tri_neg, sum_race_tri_neg)  
}

process_data(race_tri_neg_data, triple_neg_metadata_path, race_tri_neg_metadata_path, expr_dir, race_tri_neg_expr_data_path)

race_tri_neg_files <- list.files(race_tri_neg_metadata_path, full.names = TRUE)
for (file in race_tri_neg_files) {
    new_metadata <- read_tsv(file) |> 
        filter(tri_neg_status == "tri_neg")

    write_tsv(new_metadata, file.path(race_tri_neg_metadata_path, paste0(new_metadata$Dataset_ID[1], ".tsv")))
}

all_race <- race_tri_neg |> 
    filter(race %in% c("Black", "White")) |> 
    pivot_wider(names_from = tri_neg_status, values_from = n) |>
    drop_na()

big_race_final <- all_race |> 
  pivot_wider(names_from = race, values_from = c(tri_neg, non_tri_neg)) |> 
  drop_na() |>
  janitor::adorn_totals()

# calculate total number of samples
common_datasets <- list(big_race_data, big_ER_data, big_PR_data, big_HER2_data)
big_data <- reduce(common_datasets, full_join, by = "Dataset_ID") |> 
    arrange((Dataset_ID))

total_samples <- tibble(Dataset_ID = character(), num_of_samples = integer())

for (metafile in meta_dir_paths) {
    gseID <- metafile |> 
        basename() |> 
        file_path_sans_ext()
  
    if (gseID %in% big_data$Dataset_ID) {
        metadata <- read_tsv(metafile) 

        # Extract dataset ID and sample count
        Data_ID <- metadata$Dataset_ID[1]
        sample_count  <- nrow(metadata)

        # Create a tibble with correct column names
        vector_sums <- tibble(Dataset_ID = Data_ID, num_of_samples = sample_count)

        # Add to the total_samples tibble
        total_samples <- bind_rows(total_samples, vector_sums)
    }
}

# Calculate totals
total_sample_df <- total_samples |> 
  adorn_totals()

# write summary information to file
write_tsv(big_race_data, file.path(summary_data_path, "race_data.tsv"))
write_tsv(big_ER_data, file.path(summary_data_path, "ER_data.tsv"))
write_tsv(big_PR_data, file.path(summary_data_path, "PR_data.tsv"))
write_tsv(big_HER2_data, file.path(summary_data_path, "HER2_data.tsv"))
write_tsv(big_tri_neg_data, file.path(summary_data_path, "triple_negative_data.tsv"))
write_tsv(big_race_final, file.path(summary_data_path, "race_tri_neg_data.tsv"))
write_tsv(total_sample_df, file.path(summary_data_path, "Sample_numbers.tsv"))

# extra code to get number of samples for race_tri_neg
# data_dir <- "/inwosu/Meta_Analysis/Data/analysis_ready_metadata/race_tri_neg/"
# data_dir_paths <- list.files(data_dir, full.names = T)
# race_tri_neg_data <- tibble()

# for (i in seq_along(data_dir_paths)) {
  
#   sum_race <- data.frame(Dataset_ID = character(0), race = character(0), n = numeric(0))
  
#   file <- data_dir_paths[i]
#   gseID <- file |> basename() |> file_path_sans_ext()
  
#   metadata <- read_tsv(file) |>
#     mutate(across(everything(), as.character))
  
#   if ("race" %in% names(metadata)) {
#     sum_race <- metadata |> group_by(race) |> tally() |> mutate(Dataset_ID = gseID) |> relocate(Dataset_ID)
#   }
#   race_tri_neg_data <- rbind(race_tri_neg_data, sum_race)
# }

# big_race_tri_neg_data <- race_tri_neg_data |>
#   filter(race %in% c("Black", "White")) |>
#   pivot_wider(names_from = race, values_from = n) |>
#   mutate(proportion = (Black/White)) |>
#   janitor::adorn_totals()