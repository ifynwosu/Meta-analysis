
#Define input directory and list of dataset files
input_dir <- "/inwosu/Meta_Analysis/Data/prelim_expression_data/"
file_list  <- list.files(input_dir, full.names = T)[1:2]

# read in data into a list
gene_list <- list()

dataset_list <- map(file_list, function(file) {
  
  df <- read_tsv(file) |> 
    select(-c("Entrez_Gene_ID", "Ensembl_Gene_ID", "Chromosome", "Gene_Biotype")) |> 
    distinct(HGNC_Symbol, .keep_all = TRUE)
  
  Dataset_ID = df[[1, 1]]
  gene_list[[length(gene_list) + 1]] <<- df$HGNC_Symbol
  
  transposed_df <- df |> 
    t() |>
    row_to_names(2, remove_row = TRUE, remove_rows_above = TRUE) |> 
    as.data.frame() |> 
    rownames_to_column("Sample_ID") 
  
  final_df <- transposed_df |> 
    mutate(Dataset = Dataset_ID, .before = "Sample_ID")
  
  return(final_df)
})

common_genes <- tibble(Reduce(intersect, gene_list))
# get list of common genes source_col = "dataset_id"

#combine_common_cols <- function(df_list) {
combine_common_cols_existing_source <- function(df_list, source_col = "Dataset") {
  # Get common column names across all datasets
  common_cols <- reduce(dataset_list, ~ intersect(names(.x), names(.y)))
  common_cols <- union(source_col, common_cols)
  
  # Select only those common columns from each dataset
  df_list_common <- map(dataset_list, ~ select(.x, all_of(common_cols)))
  
  # Add source column and bind together
  df_combined <- bind_rows(df_list_common)
  
  return(df_combined)
}

combined_df <- combine_common_cols_existing_source(dataset_list) #|> 
  #na.omit() 

pca_data <- combined_df |> 
  select(-c(Sample_ID, Dataset)) 

pca_numeric <- map_df(pca_data, as.numeric)
pca_result <- PCA(pca_numeric, scale.unit = TRUE, graph = FALSE)
pc2 <- prcomp(pca_numeric, center=TRUE, scale=TRUE)
pc3<- as.data.frame(pc2$x)

ggplot(pc3, aes(x=PC1, y=PC2)) + 
  geom_point(aes(col=combined_df$Dataset)) +
  theme_bw(base_size=18) +
  labs(col = "Iris Species") # This changes the legend title

fviz_eig(pca_result, addlabels = TRUE, ylim = c(0, 50))

fviz_pca_biplot(pca_result,
                geom.ind = "point",
                col.ind = combined_df$Dataset,
                pointshape = 21,
                palette = "jco",
                addEllipses = TRUE,
                legend.title = "Dataset")

df1 <- tibble(x = 1:3)
df2 <- tibble(x = 3:5)
