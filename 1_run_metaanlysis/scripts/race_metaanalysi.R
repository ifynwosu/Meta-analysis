# Define the vector
my_vector <- c(1, 2, 3, 4, 5, 6, 7)

# Generate all combinations of size 6
combinations <- combn(my_vector, 6)

meta_dir <- "/inwosu/Meta_Analysis/Data/race_metadata"
meta_dir_paths <- list.files(meta_dir, full.names = T)
meta_list <- list()
#a <- as.data.frame(meta_dir_paths)

expr_dir  <- "/inwosu/Meta_Analysis/Data/race_expression_data"
expr_dir_paths <- list.files(expr_dir, full.names = T)
expr_list <- list()

for (i in seq_along(combinations[,1])) {
#for (i in seq_along(expr_dir_paths)) {
  metadata <- read_tsv(file = meta_dir_paths[combinations[,1][i]]) %>%
    column_to_rownames(var = "Sample_ID") %>%
    filter(race %in% c("Black", "White")) %>%
    filter(!is.na(race))
  meta_list[[i]] <- metadata

  expr_data <- read_tsv(file = expr_dir_paths[combinations[,1][i]]) %>%
    dplyr::select(-c("Dataset_ID","Entrez_Gene_ID", "Chromosome", "Ensembl_Gene_ID", "Gene_Biotype")) %>%
    distinct(HGNC_Symbol, .keep_all = TRUE) %>%
    column_to_rownames(var = "HGNC_Symbol") %>%
    na.omit() %>%
    as.matrix()
  expr_list[[i]] <- expr_data
}

 # account for batch effects?

phenoGroups = rep("race", (length(meta_dir_paths) - 1))
phenoCases = rep(list("Black"), (length(meta_dir_paths) - 1))
phenoControls = rep(list("White"), (length(meta_dir_paths) - 1))

meta_object <- createObjectMA(listEX = expr_list, 
               listPheno = meta_list, 
               namePheno = phenoGroups, 
               expGroups = phenoCases, 
               refGroups = phenoControls)

names(meta_object) <- str_replace_all(expr_dir_paths, c("/inwosu/Meta_Analysis/Data/race_expression_data/" = "", ".tsv.gz" = ""))

resultsMA <- metaAnalysisDE(meta_object, typeMethod = "REM")
new_results <- resultsMA %>% filter(FDR < 0.05)

heatmap_values <- draw_Heatmap(objectMA = meta_object, 
                               resMA = resultsMA,
                               typeMethod = "REM",
                               scaling = "zscor", 
                               regulation = "all",
                               numSig = 500,
                               na_col = "black",
                               legend = T,
                               case = "Black",
                               control = "White",
                               title = "Race")

Top_genes = as.data.frame(rownames(heatmap_values))

write_tsv(Top_genes, file.path("/inwosu/Meta_Analysis/Data/", "race_genes.tsv"))

effects <- calculateES(meta_object)
head(effects$ES)


# select genes based on 6 datasets
# FDR less than 0.05
# 7th dataset, narrow the dataset down to the selected genes and save that file with those gene expression values and race to tsv file
# save race values
# sklearn leave one out cross validation
# AUCROC value for the cross validation
# rotate
 
combinations[,2][2]

