library(tidyverse)

wilcox_result = tibble()

plot_graph <- function(variable, variable_plot) {
    
    results_df <- read_tsv(paste0("/Data/cross_validation_results/accuracy_file_result_", variable, ".txt"), col_names = F)
    results_df[results_df == "Meta genes"] <- "Meta-analysis genes"  
    names(results_df) <- c("Dataset_ID", "num_genes", "gene_predictor", "Balanced accuracy") 

    filtered_results_meta <- subset(results_df, num_genes == "20" & gene_predictor == "Meta-analysis genes") |>
      select(`Balanced accuracy`) |>
      unlist()
    filtered_results_random <- subset(results_df, num_genes == "20" & gene_predictor == "Random genes") |>
      select(`Balanced accuracy`) |>
      unlist()
    
    wil_result <- wilcox.test(filtered_results_meta, filtered_results_random, alternative = "greater") 
    new_vector <- tibble(
      variable = variable, 
      median_meta = median(filtered_results_meta),
      median_random = median(filtered_results_random),
      p_value = wil_result$p.value     
    )
    wilcox_result <<- bind_rows(wilcox_result, new_vector)
    
    png(paste0("/Data/cross_validation_results/", variable_plot, ".png"), width = 18, height = 10, units = "in", res = 300)

    results_df %>%
        ggplot(aes(x = factor(num_genes), y = `Balanced accuracy`, fill = factor(gene_predictor))) +
        theme(text = element_text(size=30)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) + 
        #geom_jitter(color = "black", size = 0.3, alpha = 0.4) + 
        geom_jitter(position = position_jitter(width = 0.2), size = 0.8,  alpha = 1) +        
        facet_wrap(~gene_predictor) +
        guides(fill = "none") +
        theme(axis.title.x=element_text(vjust=-3)) +
        theme(axis.title.y=element_text(vjust=3)) +
        theme(plot.margin = unit(c(1,1,1,1), "cm")) +
        xlab("Number of genes used for cross-validation predictions")
}

plot_graph("race", "race_plot")
plot_graph("ER_status", "ER_plot")
plot_graph("PR_status", "PR_plot")
plot_graph("HER2_status", "HER2_plot")
plot_graph("tri_neg_status", "tri_neg_plot")
   
# write_tsv(big_df, file.path(paste0(results_dir, variable, "_wil_result.tsv")))
write_tsv(wilcox_result, file.path("/Data/cross_validation_results/", "Wilcox_results.tsv"))
