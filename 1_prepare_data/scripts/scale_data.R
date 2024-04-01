
expression_data <- "/Data/prelim_expression_data"
output_dir <- "/Data/analysis_ready_expression_data/"

create_directory <- function(directory_path) {
    if (!dir.exists(directory_path)) {
        dir.create(directory_path, recursive = TRUE)
    }
}

create_directory(output_dir)

file_paths <- list.files(expression_data, full.names = T)

special_cases <- c("GSE62944_Tumor", "GSE62944_Normal", "ABiM.100", "ABiM.405", "Normal.66", "OSLO2EMIT0.103", "SCANB.9206")

rnaSeq_data <- function(expr_file) {
    
    #add "1" to each value in matrix
    expr_file <- sweep(expr_file, 1, 1, "+")
    
    #log2 transform values
    expr_file = log2(expr_file) 
}  

microarray_data <- function(expr_file) {
  
  #add minimum of each row to dataframe to use as a check
  expr_file$min <- apply(expr_file, 1, min)
   
  #subtract minimum value from each element in row
  expr_file <- sweep(expr_file, 1, expr_file$min ,"-")

  #add "1" to each value in matrix
  expr_file <- sweep(expr_file, 1, 1, "+") 

  expr_file <- expr_file |>
     dplyr::select(-min)
}  

for (file in file_paths) {
    cat("\n")
    print(paste0("Reading in ", file, "!"))
    cat("\n")
    
    expr_data <- read_tsv(file) %>%
        dplyr::select(-c("Dataset_ID", "Entrez_Gene_ID", "Chromosome", "HGNC_Symbol", "Gene_Biotype")) |>
        # dplyr::filter(!(if_all(everything(), ~. == 0))) |> #remove rows with all 0
        column_to_rownames(var = "Ensembl_Gene_ID") 
    
    threshold <- 0
    index <- rowSums(expr_data > threshold) 
    expr_data <- expr_data[index,]

    # expr_data <- expr_data[rowSums(expr_data[]) != 0, ]
    
    ID <- basename(file)
    ID <- gsub(".tsv.gz", "", ID)
    
    if (ID == "ICGC_KR") {
        gene_df <- expr_data 
    }   else if (ID %in% special_cases) {
            gene_df <- rnaSeq_data(expr_data) 
    }   else {
            gene_df <- microarray_data(expr_data) 
    }   

    gene_df <- gene_df |>
        rownames_to_column(var = "Ensembl_Gene_ID")

    write_tsv(gene_df, file.path(paste0(output_dir, ID, ".tsv.gz")))
    
}