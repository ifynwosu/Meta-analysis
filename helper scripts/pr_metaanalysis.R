
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 1000)

meta_dir <- "/inwosu/Meta_Analysis/Data/pr_metadata"
meta_dir_paths <- list.files(meta_dir, full.names = T)
meta_list <- list()
#a <- as.data.frame(meta_dir_paths)

expr_dir  <- "/inwosu/Meta_Analysis/Data/pr_expression_data"
expr_dir_paths <- list.files(expr_dir, full.names = T)
expr_list <- list()
#a <- as.data.frame(expr_dir_paths)

for (i in 1:3) {
#for (i in seq_along(expr_dir_paths)) {
  metadata <- read_tsv(file = meta_dir_paths[i]) %>%
    column_to_rownames(var = "Sample_ID") 
  meta_list[[i]] <- metadata
  
  expr_data <- read_tsv(file = expr_dir_paths[i]) %>%
    dplyr::select(-c("Dataset_ID","Entrez_Gene_ID", "Chromosome", "Ensembl_Gene_ID", "Gene_Biotype")) %>%
    distinct(HGNC_Symbol, .keep_all = TRUE) %>%
    column_to_rownames(var = "HGNC_Symbol") %>%
    na.omit() %>%
    as.matrix()
  expr_list[[i]] <- expr_data
}

#account for batch effects?

phenoGroups = rep("PR_Status", length(meta_dir_paths))
phenoCases = rep(list("positive"), length(meta_dir_paths))
phenoControls = rep(list("negative"), length(meta_dir_paths))

meta_object <- createObjectMA(listEX = expr_list, 
                              listPheno = meta_list, 
                              namePheno = phenoGroups, 
                              expGroups = phenoCases, 
                              refGroups = phenoControls)

names(meta_object) <- str_replace_all(expr_dir_paths, c("/inwosu/Meta_Analysis/Data/pr_expression_data/" = "", ".tsv.gz" = ""))

heterogeneityTest(meta_object)

resultsMA <- metaAnalysisDE(meta_object, typeMethod = "REM")
new_results <- resultsMA %>% filter(FDR < 0.05)

heatmap_values <- draw_Heatmap(objectMA = meta_object, 
            resMA = new_results,
            typeMethod = "REM",
            scaling = "zscor", 
            regulation = "all",
            numSig = 500,
            na_col = "black",
            legend = T,
            case = "PR_positive",
            control = "PR_negative",
            title = "PR_Status")

Top_genes = as.data.frame(rownames(heatmap_values))

write_tsv(Top_genes, file.path("/inwosu/Meta_Analysis/Data/", "PR_genes.tsv"))

