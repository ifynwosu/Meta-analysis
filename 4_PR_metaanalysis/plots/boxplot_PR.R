
accuracy_file <- read_tsv("/inwosu/Meta_Analysis/accuracy_file_ER.txt", col_names = F) 
names(accuracy_file) <- c("Dataset_ID", "meta_gene", "meta_accuracy", "random_genes", "random_accuracy", "num_genes") 

new_file <- accuracy_file |>
  mutate(new_col = ceiling(row_number() / 5)) |>
  group_by(new_col, Dataset_ID, meta_accuracy, random_accuracy, num_genes) |>
  summarise(meta = paste0(meta_gene, collapse = ","), random = paste0(random_genes, collapse = ",")) |>
  ungroup()
meta_unique <- unique(accuracy_file$meta_accuracy) |> tibble()

results_df <- read_tsv("/inwosu/Meta_Analysis/accuracy_file_ER_result.txt", col_names = F) 
names(results_df) <- c("Dataset_ID", "num_genes", "filtering_type", "accuracy_score") 


results_df %>%
  ggplot(aes(x = factor(num_genes), y = accuracy_score, fill = factor(filtering_type))) +
  geom_jitter(color = "black", size = 0.3, alpha = 0.4) +
  geom_boxplot(alpha = 0.5) + 
  facet_wrap(~filtering_type) +
  xlab("number of genes used in prediction")



