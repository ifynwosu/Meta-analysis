library(tidyverse)
library(tools)
library(janitor)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000)

# function to create directory
create_directory <- function(directory_path) {
    if (!dir.exists(directory_path)) {
        dir.create(directory_path, recursive = TRUE)
    }
    return(directory_path)
}

data_race <- create_directory("/Data/cross_validation_data/race")
data_ER <- create_directory("/Data/cross_validation_data/ER_status")
data_PR <- create_directory("/Data/cross_validation_data/PR_status")
data_HER2 <- create_directory("/Data/cross_validation_data/HER2_status")
data_trip_neg <- create_directory("/Data/cross_validation_data/tri_neg_status")
data_race_trip_neg <- create_directory("/Data/cross_validation_data/race_tri_neg")

meta <- "/Data/analysis_ready_metadata/"
expr <- "/Data/analysis_ready_expression_data/"

prepare_data <- function(path, variable, filter_variable, filter_values) {

    metadata_dir <- paste0(meta, variable)
    file_paths_meta <- list.files(metadata_dir, full.names = T)

    expr_dir <- paste0(expr, variable)
    file_paths_expr <- list.files(expr_dir, full.names = T)
    
    name_list <- list()
        
        for (i in 1:length(file_paths_expr)) {
            expr_file <- (file_paths_expr[i])
            expr_data <- read_tsv(expr_file) 
            name_list[[i]] <- expr_data$HGNC_Symbol
        }
        
    common_genes <- Reduce(intersect, name_list)

    data_tally <- list()
    
    for (i in 1:length(file_paths_meta)) {

        meta_file <- (file_paths_meta[i])
        meta_data <- read_tsv(meta_file)
        
        identifiers <- meta_data |>
            dplyr::select(Sample_ID, {{filter_variable}}) |> 
            as_tibble()
        
        filename <- meta_file |> 
            basename() |> 
            file_path_sans_ext()
        
        expr_file <- (file_paths_expr[i])
        expr_data <- read_tsv(expr_file) 

        transposed_data <- expr_data |> 
            distinct(HGNC_Symbol, .keep_all = TRUE) |>
            pivot_longer(cols= -1) |>
            pivot_wider(names_from = HGNC_Symbol, values_from = value) |>
            rename(Sample_ID = name) |>
            dplyr::select(Sample_ID, all_of(common_genes))

        joined_data <- inner_join(identifiers, transposed_data) |>
            dplyr::select(-Sample_ID) |>
            filter({{filter_variable}} %in% filter_values) |>
            filter(!is.na({{filter_variable}})) |>
            mutate(across(where(is.numeric), ~ as.numeric(scale(.x)))) # https://stackoverflow.com/questions/70377294/why-does-mutateacross-with-scale-add-1-to-the-column-header
        
        cat("\n")
        print(paste0("writing ",  filename, " to file"))
        cat("\n")

        write_tsv(joined_data, file.path(path, paste0(filename, ".tsv")))
        
        num_par <- joined_data |> 
            group_by({{filter_variable}}) |>
            tally()

        data_tally[[i]] <- num_par
        names(data_tally)[i] <- filename
    }

    summary <- bind_rows(data_tally, .id = "Dataset")
    write_tsv(summary, file.path(paste0("/Data/metaanalysis_results/", variable, "_summary.tsv")))
}

# prepare_data <- function(path, variable, filter_variable, filter_values)
prepare_data(data_race, "race", race, c("Black", "White"))
prepare_data(data_ER, "ER_status", ER_status, c("positive", "negative"))
prepare_data(data_PR, "PR_status", PR_status, c("positive", "negative"))
prepare_data(data_HER2, "HER2_status", HER2_status, c("positive", "negative"))
prepare_data(data_trip_neg, "tri_neg_status", tri_neg_status, c("tri_neg", "non_tri_neg"))
prepare_data(data_race_trip_neg, "race_tri_neg", race, c("Black", "White"))
