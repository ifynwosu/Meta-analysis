
# directory where we are going to save the results
meta_results_dir <- "/Data/results/race/"

if (!dir.exists(meta_results_dir)) {
  dir.create(meta_results_dir, recursive = TRUE)
}

# directory where metadata is stored
meta_dir <- "/Data/metadata/race"
meta_dir_paths <- list.files(meta_dir, full.names = T)
meta_list <- list()

# directory where expression data is stored
expr_dir  <- "/Data/expression_data/race"
expr_dir_paths <- list.files(expr_dir, full.names = T)
expr_list <- list()

# get common genes to use in meta-analysis
gene_list <- list()

for (i in 1:length(expr_dir_paths)) {
  expr_file <- (expr_dir_paths[i])
  expr_data <- read_tsv(expr_file) 
  gene_list[[i]] <- expr_data$HGNC_Symbol
}

common_genes <- tibble(Reduce(intersect, gene_list))
names(common_genes) <- "HGNC_Symbol"

# Define a vector to save combinations for use in generating meta-analysis list
my_vector <- c(length(meta_dir_paths)) 

# Generate all combinations of size (n - 1) because we are leaving one data set out
combinations <- combn(my_vector, (length(meta_dir_paths) - 1))

# create vectors to store values
list_all_metadata <- list()
list_all_exprdata <- list()

# prepare metadata and expression data and save into list 
for (i in 1:ncol(combinations)) {
  for (j in 1:nrow(combinations)) {
    metadata <- read_tsv(file = meta_dir_paths[combinations[,i][j]]) |>
                column_to_rownames(var = "Sample_ID") |>
                filter(race %in% c("Black", "White")) |>
                filter(!is.na(race))
    meta_list[[j]] <- metadata
    names(meta_list)[j] <- unique(meta_list[[j]][["Dataset_ID"]])
  
    expr_data <-read_tsv(file = expr_dir_paths[combinations[,i][j]]) |>
      dplyr::select(-c("Dataset_ID", "Entrez_Gene_ID", "Chromosome", "Ensembl_Gene_ID", "Gene_Biotype")) |>
      distinct(HGNC_Symbol, .keep_all = TRUE) |>
      inner_join(common_genes) |>
      column_to_rownames(var = "HGNC_Symbol") |>
      as.matrix()
    expr_list[[j]] <- expr_data
    names(expr_list)[j] <- unique(meta_list[[j]][["Dataset_ID"]])
  }
  
  list_all_metadata[[i]] <- meta_list
  list_all_exprdata[[i]] <- expr_list
}

# vector with names of all data sets
dataset_vector = str_replace_all(expr_dir_paths, c("/Data/expression_data/race/" = "", ".tsv.gz" = ""))

# create vectors for phenotype data
phenoGroups = rep("race", length(combinations[,1]))
phenoCases = rep(list("Black"), length(combinations[,1]))
phenoControls = rep(list("White"), length(combinations[,1]))

# create vector to store results from meta-analysis
meta_analysis_genes <- list()

# Perform meta-analysis
for (i in seq_along(list_all_exprdata)) {
  meta_object <- createObjectMA(listEX = list_all_exprdata[[i]], 
                              listPheno = list_all_metadata[[i]], 
                              namePheno = phenoGroups, 
                              expGroups = phenoCases, 
                              refGroups = phenoControls)
  
  meta_results <- metaAnalysisDE(meta_object, typeMethod = "REM", missAllow = 0.1, proportionData = 0.9) 
  
  new_results <- meta_results |>
    filter(FDR < 0.05) |>
    arrange(desc(abs(Com.ES))) |> 
    as_tibble(rownames = NA) |> 
    rownames_to_column("Gene")

  meta_analysis_genes[[i]] <- new_results$Gene
  
  missing_dataset <- setdiff(dataset_vector, names(meta_object))
  write_tsv(new_results, file.path(meta_results_dir, paste0("meta_results_without_", missing_dataset, ".tsv")))
}

common_meta_genes <- tibble(Reduce(intersect, meta_analysis_genes))
names(common_meta_genes) <- "Genes"

write_tsv(common_meta_genes, file.path("/Data/results/", "race_common_genes.tsv"))

# This section is for meta analysis on the full data sets.
for (i in seq_along(expr_dir_paths)) {
  metadata <- read_tsv(file = meta_dir_paths[i]) |>
    column_to_rownames(var = "Sample_ID") |>
    filter(race %in% c("Black", "White")) |>
    filter(!is.na(race))
  meta_list[[i]] <- metadata
  names(meta_list)[i] <- metadata[["Dataset_ID"]][[1]]

  expr_data <- read_tsv(file = expr_dir_paths[i]) |>
    dplyr::select(-c("Dataset_ID", "Entrez_Gene_ID", "Chromosome", "Ensembl_Gene_ID", "Gene_Biotype")) |>
    distinct(HGNC_Symbol, .keep_all = TRUE) |>
    column_to_rownames(var = "HGNC_Symbol") |>
    as.matrix()
  expr_list[[i]] <- expr_data
  names(expr_list)[i] <- metadata[["Dataset_ID"]][[1]]
}

phenoGroups = rep("race", length(meta_dir_paths))
phenoCases = rep(list("Black"), length(meta_dir_paths))
phenoControls = rep(list("White"), length(meta_dir_paths))

full_meta_object <- createObjectMA(listEX = expr_list,
                              listPheno = meta_list,
                              namePheno = phenoGroups,
                              expGroups = phenoCases,
                              refGroups = phenoControls)

full_meta_results <- metaAnalysisDE(full_meta_object, typeMethod = "REM", missAllow = 0.1, proportionData = 0.9)

new_full_results <- full_meta_results |>
    filter(FDR < 0.05) |>
    arrange(desc(abs(Com.ES))) |>
    as_tibble(rownames = NA) |>
    rownames_to_column("Gene")

write_tsv(new_full_results, file.path("/Data/results/", "race_full_meta_results.tsv"))

source("/prepare_data/functions/draw_Heatmap.R")
png("/Data/results/race_heatmap.png", width = 24, height = 12, units = "in", res = 300)
heatmap_values <- draw_Heatmap(objectMA = full_meta_object,
                               resMA = full_meta_results,
                               typeMethod = "REM",
                               scaling = "zscor",
                               fdrSig = 0.05,
                               regulation = "all",
                               numSig = 1000,
                               case = "Black",
                               control = "White",
                               title = "Race")
dev.off()

top_genes = as.data.frame(rownames(heatmap_values)) #check if top_genes is the same as new_result 
write_tsv(top_genes, file.path("/Data/results/", "race_heatmap_genes.tsv"))

# effects <- calculateES(meta_object)
# head(effects$ES)