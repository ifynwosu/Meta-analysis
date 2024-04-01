library(tidyverse)

plot_graph <- function(variable, variable_plot) {
    
    results_df <- read_tsv(paste0("/Data/cross_validation_results/accuracy_file_result_", variable, ".txt"), col_names = F) 
    names(results_df) <- c("Dataset_ID", "num_genes", "gene_predictor", "Balanced accuracy") 

    png(paste0("/Data/cross_validation_results/", variable_plot, ".png"), width = 18, height = 10, units = "in", res = 300)

    results_df %>%
        ggplot(aes(x = factor(num_genes), y = `Balanced accuracy`, fill = factor(gene_predictor))) +
        theme(text = element_text(size=30)) +
        #geom_jitter(color = "black", size = 0.3, alpha = 0.4) + 
        geom_jitter(position = position_jitter(width = 0.2), size = 0.3,  alpha = 0.5) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
        facet_wrap(~gene_predictor) +
        guides(fill = "none") +
        xlab("number of genes used in prediction")

}

plot_graph("race", "race_plot")
plot_graph("ER_status", "ER_plot")
plot_graph("PR_status", "PR_plot")
plot_graph("HER2_status", "HER2_plot")
plot_graph("tri_neg_status", "tri_neg_plot")
   


# ggsave("/inwosu/Meta_Analysis/Data/cross_validation_results/PR_plot.png")
# 
# accuracy_file <- read_tsv("/inwosu/Meta_Analysis/accuracy_file_ER.txt", col_names = F) 
# names(accuracy_file) <- c("Dataset_ID", "meta_gene", "meta_accuracy", "random_genes", "random_accuracy", "num_genes") 
# 
# new_file <- accuracy_file |>
#   mutate(new_col = ceiling(row_number() / 5)) |>
#   group_by(new_col, Dataset_ID, meta_accuracy, random_accuracy, num_genes) |>
#   summarise(meta = paste0(meta_gene, collapse = ","), random = paste0(random_genes, collapse = ",")) |>
#   ungroup()
# meta_unique <- unique(accuracy_file$meta_accuracy) |> tibble()