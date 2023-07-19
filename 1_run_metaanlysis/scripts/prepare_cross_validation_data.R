
cross_validation_data_dir <- "/inwosu/Meta_Analysis/Data/cross_validation"

if (!dir.exists(cross_validation_data_dir)) {
  dir.create(cross_validation_data_dir)
}

metadata_dir <- "/inwosu/Meta_Analysis/Data/race_metadata"
file_paths_meta <- list.files(metadata_dir, full.names = T)

expr_dir <- "/inwosu/Meta_Analysis/Data/race_expression_data"
file_paths_expr <- list.files(expr_dir, full.names = T)

# get common genes
name_list <- list()

for (i in 1:length(file_paths_meta)) {
  expr_file <- (file_paths_expr[i])
  expr_data <- read_tsv(expr_file) 
  name_list[[i]] <- expr_data$HGNC_Symbol
}

common_genes <- Reduce(intersect, name_list)

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
    
  transposed_data <- expr_data |> 
    dplyr::select(-c("Dataset_ID", "Entrez_Gene_ID", "Chromosome", "Ensembl_Gene_ID", "Gene_Biotype")) |>
    distinct(HGNC_Symbol, .keep_all = TRUE) |>
    t() |>
    as_tibble(rownames = NA) |> 
    rownames_to_column()
  
  renamed_data <- row_to_names(transposed_data, 1, remove_row = TRUE, remove_rows_above = TRUE) |>
    rename(Sample_ID = HGNC_Symbol) |>
    select(Sample_ID, all_of(common_genes))
   
  joined_data <- inner_join(identifiers, renamed_data) |>
    select(- Sample_ID) |>
    filter(race %in% c("Black", "White")) |>
    filter(!is.na(race))
  
  write_tsv(joined_data, file.path(cross_validation_data_dir, paste0(filename, ".tsv")))
}

  

