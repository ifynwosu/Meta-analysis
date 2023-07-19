
cross_validation_results_dir <- "/inwosu/Meta_Analysis/Data/cross_validation"

if (!dir.exists(cross_validation_results_dir)) {
  dir.create(cross_validation_results_dir)
}

meta_dir <- "/inwosu/Meta_Analysis/Data/race_metadata"
# meta_dir <- "/inwosu/Meta_Analysis/Data/meta_results"
file_paths_meta <- list.files(meta_dir, full.names = T)

expr_dir <- "/inwosu/Meta_Analysis/Data/race_expression_data"
file_paths_expr <- list.files(expr_dir, full.names = T)

results_dir <- "/inwosu/Meta_Analysis/Data/meta_results"
file_paths_results <- list.files(results_dir, full.names = T)

for (i in 1:length(file_paths_meta)) {
  
  meta_file <- (file_paths_meta[i])
  meta_data <- read_tsv(meta_file)
  
  identifiers <- meta_data %>%
    select(Sample_ID, race) %>% as_tibble()
  
  filename <- meta_file %>% basename() %>% file_path_sans_ext()
  
  expr_file <- (file_paths_expr[i])
  expr_data <- read_tsv(expr_file) 
  
  result_file <- (file_paths_results[i])
  result_data <- read_tsv(result_file)
  relevant_genes <- result_data %>%
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
    select(- Sample_ID) |>
    filter(race %in% c("Black", "White")) |>
    filter(!is.na(race))
  
  write_csv(joined_data, file.path(cross_validation_results_dir, paste0(filename, ".csv")))
}

