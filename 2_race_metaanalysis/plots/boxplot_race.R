
accuracy_file <- read_tsv("/inwosu/Meta_Analysis/accuracy_file_race.txt", col_names = F) 
names(accuracy_file) <- c("Dataset_ID", "meta_gene", "meta_accuracy", "random_genes", "random_accuracy", "num_genes") 

new_file <- accuracy_file |>
  mutate(new_col = ceiling(row_number() / 5)) |>
  group_by(new_col, Dataset_ID, meta_accuracy, random_accuracy, num_genes) |>
  summarise(meta = paste0(meta_gene, collapse = ","), random = paste0(random_genes, collapse = ",")) |>
  ungroup()
meta_unique <- unique(accuracy_file$meta_accuracy) |> tibble()

race_results <- read_tsv("/inwosu/Meta_Analysis/accuracy_file_race_result.txt", col_names = F) 
names(race_results) <- c("Dataset_ID", "num_genes", "gene_predictor", "Balanced accuracy") 


race_results %>%
  ggplot(aes(x = factor(num_genes), y = `Balanced accuracy`, fill = factor(gene_predictor))) +
  geom_jitter(color = "black", size = 0.3, alpha = 0.4) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  facet_wrap(~gene_predictor) +
  guides(fill = "none", ) +
  xlab("Number of genes used in prediction")

ggsave("raceplot.pdf")


