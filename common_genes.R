
file_dir  <- "/inwosu/Meta_Analysis/Data/HER2_meta_results"
file_dir_paths <- list.files(file_dir, full.names = T)

gene_list <- list()

for (i in 1:length(file_dir_paths)) {
#for (i in 1:6) {
  meta_file <- (file_dir_paths[i])
  meta_result <- read_tsv(meta_file) 
  gene_list[[i]] <- meta_result$Gene |>
    head(n= 100) # this is the number of genes that performed best in CV (will be diff for each variable)
  
}

common_genes <- tibble(Reduce(intersect, gene_list))
names(common_genes) <- "HGNC_Symbol"

# top_genes = as.data.frame(rownames(heatmap_values)) #check if top_genes is the same as new_result 
write_tsv(common_genes, file.path("/inwosu/Meta_Analysis/Data/", "common_race_genes.tsv"))

new_results <- meta_results |>
  filter(FDR < 0.05) |>
  arrange(desc(abs(Com.ES))) |> 
  head(n = 500) |>
  as_tibble(rownames = NA) |> 
  rownames_to_column("Gene")