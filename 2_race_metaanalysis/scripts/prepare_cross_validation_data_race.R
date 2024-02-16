
cross_validation_data <- "/Data/cross_validation_data/race"

if (!dir.exists(cross_validation_data)) {
  dir.create(cross_validation_data, recursive = TRUE)
}

metadata_dir <- "/Data/metadata/race/"
file_paths_meta <- list.files(metadata_dir, full.names = T)

expr_dir <- "/Data/expression_data/race/"
file_paths_expr <- list.files(expr_dir, full.names = T)

# results_dir <- "/Data/results/race/"
# file_paths_results <- list.files(results_dir, full.names = T)

# get common genes
name_list <- list()
for (i in 1:length(file_paths_expr)) {
  expr_file <- (file_paths_expr[i])
  expr_data <- read_tsv(expr_file) 
  name_list[[i]] <- expr_data$HGNC_Symbol
}
common_genes <- Reduce(intersect, name_list)

race_tally <- list()

# transpose expression file and select common genes
for (i in 1:length(file_paths_meta)) {

    meta_file <- (file_paths_meta[i])
    meta_data <- read_tsv(meta_file)
    
    identifiers <- meta_data |>
        select(Sample_ID, race) |> 
        as_tibble()
    
    filename <- meta_file |> 
        basename() |> 
        file_path_sans_ext()
    
    expr_file <- (file_paths_expr[i])
    expr_data <- read_tsv(expr_file) 

    # result_file <- (file_paths_results[i])
    # result_data <- read_tsv(result_file)
    # relevant_genes <- result_data |>
    #     pull(Gene)

    transposed_data <- expr_data |> 
        dplyr::select(-c("Dataset_ID", "Entrez_Gene_ID", "Chromosome", "Ensembl_Gene_ID", "Gene_Biotype")) |>
        distinct(HGNC_Symbol, .keep_all = TRUE) |>
        # filter(HGNC_Symbol %in% relevant_genes) |>
        t() |>
        as_tibble(rownames = NA) |> 
        rownames_to_column()

    renamed_data <- row_to_names(transposed_data, 1, remove_row = TRUE, remove_rows_above = TRUE) |>
        rename(Sample_ID = HGNC_Symbol) |>
    dplyr::select(Sample_ID, all_of(common_genes))

    joined_data <- inner_join(identifiers, renamed_data) |>
        select(-Sample_ID) |>
        filter(race %in% c("Black", "White")) |>
        filter(!is.na(race)) |>
        mutate(across(where(is.numeric), scale))
    
    cat("\n")
    print(paste0("writing ",  filename, " to file"))
    cat("\n")

    write_tsv(joined_data, file.path(cross_validation_data, paste0(filename, ".tsv")))
    
    num_par <- joined_data |> 
        group_by(race) |>
        tally()

    race_tally[[i]] <- num_par
    names(race_tally)[i] <- filename
}

race_summary <- bind_rows(race_tally, .id = "Dataset")
write_tsv(race_summary, file.path("/Data/results", "race_summary.tsv"))
