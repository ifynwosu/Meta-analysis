
new_data <- expr_data |> 
  dplyr::select(-c("Dataset_ID", "Entrez_Gene_ID", "Chromosome", "Ensembl_Gene_ID", "Gene_Biotype")) |>
  pivot_longer(cols = starts_with("MB"))


ggplot(new_data, aes(x = as.factor(HGNC_Symbol), y = value)) +
  geom_boxplot() +
  scale_fill_manual(values = c("#09E359", "#E31009")) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
