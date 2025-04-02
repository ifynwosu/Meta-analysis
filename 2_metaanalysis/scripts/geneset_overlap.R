
gene_set_name <- function(my_file, num_rows) {
  data_file <- my_file |>
    read_tsv() |>
    row_to_names(6) |>
    clean_names() |>
    select("gene_set_name") |>
    slice(1:num_rows)
  
  return(data_file)
}

ER_data <- gene_set_name("/Data/msigdb/Supplementary File 9.xlsx - ER_overlap_20.tsv", 20)
PR_data <- gene_set_name("/Data/msigdb/Supplementary File 10.xlsx - PR_overlap_20.tsv", 20)
HER2_data <- gene_set_name("/Data/msigdb/Supplementary File 11.xlsx - HER2_overlap_20_genes.tsv", 20)
TNBC_data <- gene_set_name("/Data/msigdb/Supplementary File 12.xlsx - TNBC_overlap_20_genes.tsv", 20)

gene_set_list <- list(ER_data, PR_data, HER2_data, TNBC_data)

common_gene_set <- reduce(gene_set_list, inner_join)