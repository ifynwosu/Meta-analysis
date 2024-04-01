

# create relevant directories

create_directory <- function(directory_path) {
    if (!dir.exists(directory_path)) {
        dir.create(directory_path, recursive = TRUE)
    }
    return(directory_path)
}

race_metadata_path <- create_directory("/Data/metadata/race/")
race_expression_data_path <- create_directory("/Data/expression_data/race/")

ER_metadata_path <- create_directory("/Data/metadata/ER_status/")
ER_expression_data_path <- create_directory("/Data/expression_data/ER_status/")

PR_metadata_path <- create_directory("/Data/metadata/PR_status/")
PR_expression_data_path <- create_directory("/Data/expression_data/PR_status/")

HER2_metadata_path <- create_directory("/Data/metadata/HER2_status/")
HER2_expression_data_path <- create_directory("/Data/expression_data/HER2_status/")

triple_neg_metadata_path <- create_directory("/Data/metadata/tri_neg_status/")
triple_neg_expr_data_path <- create_directory("/Data/expression_data/tri_neg_status/")


# directory where expression data is stored
expr_dir <- "/Data/prelim_expression_data/"

# directory where clean metadata will be stored
metadata_dir <- "/Data/analysis_ready_metadata/"
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
    
    file <- meta_dir_paths[i]
    gseID <- file |> basename() |> file_path_sans_ext()
    
    metadata <- read_tsv(file) %>%
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
    pivot_wider(names_from = race, values_from = n) |>
    mutate(proportion = (Black/White)) |>
    arrange(proportion) |>
    filter(between(proportion, 0.05, 20)) |>
    filter(Black >= 10, White >= 10) |>
    filter(!Dataset_ID == "GSE62944_Normal") %>%
    janitor::adorn_totals()

big_ER_data <- er_data |>
    pivot_wider(names_from = ER_status, values_from = n) |>
    mutate(proportion = (negative/positive)) |>
    arrange(proportion) |>
    filter(between(proportion, 0.05, 20)) |>
    filter(negative >= 10, positive >= 10) |>
    filter(!Dataset_ID == "GSE62944_Normal") %>%
    janitor::adorn_totals()

big_PR_data <- pr_data |> 
    pivot_wider(names_from = PR_status, values_from = n) |>
    mutate(proportion = (negative/positive)) |>
    arrange(proportion) |>
    filter(between(proportion, 0.05, 20)) |>
    filter(negative >= 10, positive >= 10) |>
    filter(!Dataset_ID == "GSE62944_Normal") %>%
    janitor::adorn_totals()

big_HER2_data <- her2_data |>
    pivot_wider(names_from = HER2_status, values_from = n) |>
    mutate(proportion = (negative/positive)) |>
    arrange(proportion) |>
    #mutate(across(proportion, round, 5)) |>
    filter(between(proportion, 0.05, 20)) |>
    filter(negative >= 10, positive >= 10) |>
    filter(!Dataset_ID == "GSE62944_Normal") %>%
    janitor::adorn_totals()

write_tsv(big_race_data, file.path("/Data/", "race_data.tsv"))
write_tsv(big_ER_data, file.path("/Data/", "ER_data.tsv"))
write_tsv(big_PR_data, file.path("/Data/", "PR_data.tsv"))
write_tsv(big_HER2_data, file.path("/Data/", "HER2_data.tsv"))

## The code block below identifies datasets with all 3 HR statuses
ER <- big_ER_data$Dataset_ID |> as_tibble()
PR <- big_PR_data$Dataset_ID |> as_tibble()
HER2 <- big_HER2_data$Dataset_ID |> as_tibble()

combined_data <- intersect(ER, PR)
combined_data <- intersect(combined_data, HER2)
names(combined_data) <- "Dataset_ID"


# copy necessary files to appropriate folders for analysis in later steps
copy_data <- function(names_file, data_from, data_to) {
    for (i in names_file) {
        dataPath <- paste0(data_from, i)
        file.copy(dataPath, data_to)
    }
}

process_data <- function(data, metadata_dir, metadata_path, expr_dir, expression_data_path) {    
    data$Dataset_ID <- paste0(data$Dataset_ID, ".tsv") # Add ".tsv" extension to dataset id    
    names_file <- data$Dataset_ID |>        
        as_tibble() |>
        dplyr::filter(row_number() <= n() - 1)          # remove last row (contains "Totals" from summing)

    copy_data(names_file, metadata_dir, metadata_path)
    
    data$Dataset_ID <- paste0(data$Dataset_ID, ".gz") # Add ".gz" extension to dataset id + .tsv from above
    names_file <- data$Dataset_ID |>
        as_tibble() |>
        dplyr::filter(row_number() <= n() - 1)

    copy_data(names_file, expr_dir, expression_data_path)
}

process_data(big_race_data, metadata_dir, race_metadata_path, expr_dir, race_expression_data_path)
process_data(big_ER_data, metadata_dir, ER_metadata_path, expr_dir, ER_expression_data_path)
process_data(big_PR_data, metadata_dir, PR_metadata_path, expr_dir, PR_expression_data_path)
process_data(big_HER2_data, metadata_dir, HER2_metadata_path, expr_dir, HER2_expression_data_path)
process_data(combined_data, metadata_dir, triple_neg_metadata_path, expr_dir, triple_neg_expr_data_path)


## identify triple negative samples
triple_negative_path <- list.files(triple_neg_metadata_path, full.names = TRUE)
trip_neg_data <- tibble()

identify_trip_neg <- function(data) {
    data %>%
        mutate(tri_neg_status = case_when(
                                          ER_status == "positive" ~ "non_tri_neg",
                                          PR_status == "positive" ~ "non_tri_neg",
                                          HER2_status == "positive" ~ "non_tri_neg",
                                          ER_status == "negative" & PR_status == "negative" & HER2_status == "negative" ~ "tri_neg")) |>
        select(Dataset_ID, Sample_ID, ER_status, PR_status, HER2_status, tri_neg_status, everything())
}



for (i in seq_along(triple_negative_path)) {    
    sum_trip_neg <- data.frame(Dataset_ID = character(0), tri_neg_status = character(0), n = numeric(0))
    
    meta_file <- (triple_negative_path[i])
    gseID <- meta_file |> basename() |> file_path_sans_ext()

    meta_data <- read_tsv(meta_file) |>
        identify_trip_neg() |>
        filter(!is.na(tri_neg_status))

    sum_trip_neg <- meta_data |> group_by(tri_neg_status) |> tally() |> mutate(Dataset_ID = gseID) |> relocate(Dataset_ID)    
    trip_neg_data <- rbind(trip_neg_data, sum_trip_neg)

    write_tsv(meta_data, file.path(triple_neg_metadata_path,  paste0(meta_data$Dataset_ID[1], ".tsv")))
}

big_tri_neg_data <- trip_neg_data |>
    pivot_wider(names_from = tri_neg_status, values_from = n) |>
    janitor::adorn_totals()

write_tsv(big_tri_neg_data, file.path("/Data/", "triple_negative_data.tsv"))

# #check data labels from beginning to the predicted data (GEO to predicted df)
# mydata <- read_tsv("/inwosu/Meta_Analysis/Data/metadata/triple_negative/GSE96058_HiSeq.tsv")
# 
# meta_data <- mydata |>
#   #replace_trip_neg() |>
#   filter(!is.na(result))
# 
# is.na(mydata$ER_status)
