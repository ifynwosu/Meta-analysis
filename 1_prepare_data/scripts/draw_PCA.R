
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000)

# Define input directory and list of dataset files
input_dir <- "/Data/temp_expression_data/"
file_list  <- list.files(input_dir, full.names = T)

# Dataset_Ids file
select_data <- read_tsv("/Data/summary/Sample_numbers.tsv")

plot_dir <- "/Data/"

PCA_list <- list()
i = 1

# create list of relevant datasets for PCA
for (file in file_list) {
  gseID <- file |>   
    basename() |> 
    file_path_sans_ext() |> 
    file_path_sans_ext()
    
  if (gseID %in% select_data$Dataset_ID) {
    PCA_list[[i]] <- file
    i=i+1
  }
}

# Process datasets into a list and transpose. Collect gene symbols too.
expr_data <- map(PCA_list, function(file) {
  df <- read_tsv(file) 
  
  gseID <- file |>   
    basename() |> 
    file_path_sans_ext() |> 
    file_path_sans_ext()
  
  hgnc_symbol <- pull(df, HGNC_Symbol)
  Dataset_ID <- gseID
 
  transposed_df <- df |>
    t() |>
    row_to_names(1) |>
    as.data.frame() |>
    rownames_to_column("Sample_ID")
  
  final_df <- transposed_df |>
    mutate(Dataset = Dataset_ID, .before = "Sample_ID")

  list(hgnc = hgnc_symbol, data = final_df)
})

# Split out HGNC lists and processed datasets
hgnc_list <- map(expr_data, "hgnc")
dataset_list <- map(expr_data, "data")

# Find common HGNC_Symbols across all datasets
common_genes <- reduce(hgnc_list, intersect)

combine_common_cols <- function(df_list, source_col = "Dataset") {
  
  common_cols <- union(source_col, common_genes)
  # Select only those common columns from each dataset
  df_list_common <- map(df_list, ~ select(.x, all_of(common_cols)))
  
  # Add source column and bind together
  df_combined <- bind_rows(df_list_common)
  
  return(df_combined)
}

combined_df <- combine_common_cols(dataset_list)

pca_data <- combined_df |> select(-Dataset)
pca_numeric <- map_df(pca_data, as.numeric)

pca_default <- prcomp(pca_numeric, center=TRUE, scale=TRUE)
pcdata<- as.data.frame(pca_default$x)

dataplot <- pcdata |> 
    ggplot(aes(x=PC1, y=PC2)) +
    geom_point(aes(col=combined_df$Dataset)) +
    theme_bw(base_size=10) +
    labs(col = "Datasets") + 
    theme(legend.key.size = unit(0.2, 'in')) + 
    guides(colour = guide_legend(ncol = 2)) +
    xlab("Plot for datasets used in the metaanalysis") 
    # theme(legend.title = element_text(size = 10), 
    #       legend.text = element_text(size = 6))

ggsave("/Data/PCA_plot.png", dataplot, width = 18, height = 10, units = "in", dpi = 300)


# pca_result <- PCA(pca_numeric, scale.unit = TRUE, graph = FALSE)
# fviz_pca_ind(pca_result,
#              geom.ind = "point",
#              col.ind = combined_df$Dataset,
#              pointshape = 21,
#              palette = "jco",
#              #addEllipses = TRUE,
#              legend.title = "Datasets")

