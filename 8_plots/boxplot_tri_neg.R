

results_df <- read_tsv("/inwosu/Meta_Analysis/Data/cross_validation_results/accuracy_file_result_tri_neg.txt", col_names = F) 
names(results_df) <- c("Dataset_ID", "num_genes", "gene_predictor", "Balanced accuracy") 

png("/inwosu/Meta_Analysis/Data/cross_validation_results/tri_neg_plot.png", width = 18, height = 10, units = "in", res = 300)
results_df %>%
  ggplot(aes(x = factor(num_genes), y = `Balanced accuracy`, fill = factor(gene_predictor))) +
  theme(text = element_text(size = 30)) +
  geom_jitter(color = "black", size = 0.3, alpha = 0.4) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
  facet_wrap(~gene_predictor) +
  guides(fill = "none") +
  xlab("number of genes used in prediction")
dev.off()









# 
# ggsave("HER2_plot.pdf")
# 
# accuracy_file <- read_tsv("/inwosu/Meta_Analysis/accuracy_file_HER2.txt", col_names = F) 
# names(accuracy_file) <- c("Dataset_ID", "meta_gene", "meta_accuracy", "random_genes", "random_accuracy", "num_genes") 
# 
# new_file <- accuracy_file |>
#   mutate(new_col = ceiling(row_number() / 5)) |>
#   group_by(new_col, Dataset_ID, meta_accuracy, random_accuracy, num_genes) |>
#   summarise(meta = paste0(meta_gene, collapse = ","), random = paste0(random_genes, collapse = ",")) |>
#   ungroup()
# meta_unique <- unique(accuracy_file$meta_accuracy) |> tibble()
# 
