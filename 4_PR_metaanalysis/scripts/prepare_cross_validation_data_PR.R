
cross_validation_data <- "/Data/cross_validation_data/PR_status"

if (!dir.exists(cross_validation_data)) {
  dir.create(cross_validation_data, recursive = TRUE)
}

metadata_dir <- "/Data/metadata/PR_status/"
file_paths_meta <- list.files(metadata_dir, full.names = T)

expr_dir <- "/Data/expression_data/PR_status/"
file_paths_expr <- list.files(expr_dir, full.names = T)

results_dir <- "/Data/results/PR_status/"
file_paths_results <- list.files(results_dir, full.names = T)

PR_status_tally <- list()

# transpose expression file and select common genes
for (i in 1:length(file_paths_meta)) {
  
    meta_file <- (file_paths_meta[i])
    meta_data <- read_tsv(meta_file)
    
    identifiers <- meta_data |>
        dplyr::select(Sample_ID, PR_status) |> 
        as_tibble()
    
    filename <- meta_file |> 
        basename() |> 
        file_path_sans_ext()
    
    expr_file <- (file_paths_expr[i])
    expr_data <- read_tsv(expr_file) 
    
    result_file <- (file_paths_results[i])
    result_data <- read_tsv(result_file)
    relevant_genes <- result_data |>
        pull(Gene)
    
    transposed_data <- expr_data |> 
        dplyr::select(-c("Dataset_ID", "Entrez_Gene_ID", "Chromosome", "Ensembl_Gene_ID", "Gene_Biotype")) |>
        distinct(HGNC_Symbol, .keep_all = TRUE) |>
        filter(HGNC_Symbol %in% relevant_genes) |>
        t() |>
        as_tibble(rownames = NA) |> 
        rownames_to_column()
    
    renamed_data <- row_to_names(transposed_data, 1, remove_row = TRUE, remove_rows_above = TRUE) |>
        rename(Sample_ID = HGNC_Symbol)
    
    joined_data <- inner_join(identifiers, renamed_data) |>
        dplyr::select(- Sample_ID) |>
        filter(PR_status %in% c("positive", "negative")) |>
        filter(!is.na(PR_status)) |>
        mutate(across(where(is.numeric), scale))

    cat("\n")
    print(paste0("writing ",  filename, " to file"))
    cat("\n")

    write_tsv(joined_data, file.path(cross_validation_data, paste0(filename, ".tsv")))
    
    num_par <- joined_data |> 
        group_by(PR_status) |>
        tally()
    
    PR_status_tally[[i]] <- num_par
    names(PR_status_tally)[i] <- filename
}

PR_summary <- bind_rows(PR_status_tally, .id = "Dataset")
write_tsv(PR_summary, file.path("/Data/results", "PR_summary.tsv"))




































# import os, sys

# # Data Processing
# import pandas as pd
# import numpy as np

# from numpy import mean
# from numpy import std
# from sklearn.model_selection import LeaveOneOut
# from sklearn.model_selection import cross_val_score
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.linear_model import LogisticRegression
# from sklearn.model_selection import KFold
# from sklearn.model_selection import StratifiedKFold

# # set directory
# script_path = "/home"
# os.chdir(script_path)

# # number of genes to use in cross validation
# n = 100

# # define methods
# # kfold = KFold(n_splits = 3)
# kfold = StratifiedKFold(n_splits = 5) # test diff values 
# cv_method = kfold
# cross_validation_data <- "/Data/cross_validation_meta"

# if (!dir.exists(cross_validation_data)) {
#   dir.create(cross_validation_data)
# }

# metadata_dir <- "/Data/race_metadata"
# file_paths_meta <- list.files(metadata_dir, full.names = T)

# expr_dir <- "/Data/race_expression_data"
# file_paths_expr <- list.files(expr_dir, full.names = T)

# results_dir <- "/Data/meta_results"
# file_paths_results <- list.files(results_dir, full.names = T)

# for (i in 1:length(file_paths_meta)) {
  
#   meta_file <- (file_paths_meta[i])
#   meta_data <- read_tsv(meta_file)
  
#   identifiers <- meta_data |>
#     select(Sample_ID, race) |> 
#     as_tibble()
  
#   filename <- meta_file |> 
#     basename() |> 
#     file_path_sans_ext()
  
#   expr_file <- (file_paths_expr[i])
#   expr_data <- read_tsv(expr_file) 
  
#   result_file <- (file_paths_results[i])
#   result_data <- read_tsv(result_file)
#   relevant_genes <- result_data |>
#     pull(Gene)
  
#   transposed_data <- expr_data |> 
#     dplyr::select(-c("Dataset_ID", "Entrez_Gene_ID", "Chromosome", "Ensembl_Gene_ID", "Gene_Biotype")) |>
#     distinct(HGNC_Symbol, .keep_all = TRUE) |>
#     filter(HGNC_Symbol %in% relevant_genes) |>
#     t() |>
#     as_tibble(rownames = NA) |> 
#     rownames_to_column()
  
#   renamed_data <- row_to_names(transposed_data, 1, remove_row = TRUE, remove_rows_above = TRUE) |>
#     rename(Sample_ID = HGNC_Symbol) 
  
#   joined_data <- inner_join(identifiers, renamed_data) |>
#     select(- Sample_ID) |>
#     filter(race %in% c("Black", "White")) |>
#     filter(!is.na(race))
  
#   write_csv(joined_data, file.path(cross_validation_data, paste0(filename, ".csv")))
# }


# # cv_method = LeaveOneOut()

# # create model
# # try logistic regression use default parameters for now
# logisticRegr = LogisticRegression()

# accuracy_file = open("accuracy_file.txt", "w")

# meta_results_dir = "Data/race_meta_results"
# cross_val_data_dir = "Data/cross_validation"

# logisticRegr.fit(x_train, y_train)


# Meta_Analysis/1_run_metaanlysis/scripts/prepare_cross_validation_data_from_meta_results.R