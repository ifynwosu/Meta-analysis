create_directory <- function(directory_path) {
  if (!dir.exists(directory_path)) {
    dir.create(directory_path, recursive = TRUE)
  }
  return(directory_path)
}

process_data <- function(variable) { 
  
  metadata_dir <- paste0(metadata_path, variable)
  file_paths_meta <- list.files(metadata_dir, full.names = T)
  
  expr_dir <- paste0(expr_data_path, variable)
  file_paths_expr <- list.files(expr_dir, full.names = T)
  
  metadata_final <- create_directory(paste0(metadata_output_path, variable))
  expr_data_final <- create_directory(paste0(expr_data_output_path, variable))
  
  for (i in seq_along(file_paths_meta)) {
      
      meta_file <- (file_paths_meta[i])
      meta_data <- read_tsv(meta_file)
       
      expr_file <- (file_paths_expr[i])
      expr_data <- read_tsv(expr_file)
       
      filename_meta <- meta_file %>% basename() %>% file_path_sans_ext()
      filename_expr <- expr_file %>% basename() %>% file_path_sans_ext()
      
      meta_path <- file.path(metadata_final, paste0(filename_meta, ".tsv"))
      expr_path <- file.path(expr_data_final, paste0(filename_expr, ".gz"))
          
      Sample_ID <- names(expr_data)[2:ncol(expr_data)]
      all_samples <- meta_data |>
        pull(Sample_ID)
      keep_samples <- intersect(Sample_ID, all_samples)
      
      clean_meta_data <- meta_data[meta_data$Sample_ID %in% keep_samples, ]
      clean_expr_data <- dplyr::select(expr_data, HGNC_Symbol, all_of(keep_samples))
      
      write_tsv(clean_meta_data, meta_path)
      write_tsv(clean_expr_data, expr_path)

  }
}

metadata_path <- "/Data/temp2_metadata/"
expr_data_path <- "/Data/temp2_expr_data/"

metadata_output_path <- "/Data/analysis_ready_metadata/"
expr_data_output_path <- "/Data/analysis_ready_expression_data/"

process_data("race")
process_data("ER_status")
process_data("PR_status")
process_data("HER2_status")
process_data("tri_neg_status")
process_data("race_tri_neg")
