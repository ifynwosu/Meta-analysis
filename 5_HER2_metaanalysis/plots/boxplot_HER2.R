
accuracy_file <- read_tsv("/inwosu/Meta_Analysis/accuracy_file_HER2.txt", col_names = F) 
names(accuracy_file) <- c("Dataset_ID", "meta_gene", "meta_accuracy", "random_genes", "random_accuracy", "num_genes") 

new_file <- accuracy_file |>
  mutate(new_col = ceiling(row_number() / 5)) |>
  group_by(new_col, Dataset_ID, meta_accuracy, random_accuracy, num_genes) |>
  summarise(meta = paste0(meta_gene, collapse = ","), random = paste0(random_genes, collapse = ",")) |>
  ungroup()
meta_unique <- unique(accuracy_file$meta_accuracy) |> tibble()

#start to think about 4 to 6 figures

results_df <- read_tsv("/inwosu/Meta_Analysis/accuracy_file_HER2_result.txt", col_names = F) 
names(results_df) <- c("Dataset_ID", "num_genes", "gene_predictor", "accuracy score") 

results_df %>%
  ggplot(aes(x = factor(num_genes), y = `accuracy score`, fill = factor(gene_predictor))) +
  geom_jitter(color = "black", size = 0.3, alpha = 0.4) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  facet_wrap(~gene_predictor) +
  guides(fill = "none", ) +
  xlab("number of genes used in prediction")

ggsave("HER2_plot.pdf")


