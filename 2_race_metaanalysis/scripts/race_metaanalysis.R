
# do we account for batch effects?

# directory where we are going to save the results
meta_results_dir <- "/inwosu/Meta_Analysis/Data/race_meta_results"

if (!dir.exists(meta_results_dir)) {
  dir.create(meta_results_dir, recursive = TRUE)
}

# directory where metadata is stored
meta_dir <- "/inwosu/Meta_Analysis/Data/race_metadata"
meta_dir_paths <- list.files(meta_dir, full.names = T)
meta_list <- list()

# directory where expression data is stored
expr_dir  <- "/inwosu/Meta_Analysis/Data/race_expression_data"
expr_dir_paths <- list.files(expr_dir, full.names = T)
expr_list <- list()

# get common genes
name_list <- list()

for (i in 1:length(expr_dir_paths)) {
  expr_file <- (expr_dir_paths[i])
  expr_data <- read_tsv(expr_file) 
  name_list[[i]] <- expr_data$HGNC_Symbol
}

common_genes <- tibble(Reduce(intersect, name_list))
names(common_genes) <- "HGNC_Symbol"

# Define a vector to save combinations for use in generating meta-analysis list
my_vector <- c(length(meta_dir_paths)) 

# Generate all combinations of size (n - 1) because we are leaving one data set out
combinations <- combn(my_vector, (length(meta_dir_paths) - 1))

# create vectors for phenotype data
phenoGroups = rep("race", length(combinations[,1]))
phenoCases = rep(list("Black"), length(combinations[,1]))
phenoControls = rep(list("White"), length(combinations[,1]))

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
      na.omit() |>
      as.matrix()
    expr_list[[j]] <- expr_data
  }
  
  list_all_metadata[[i]] <- meta_list
  list_all_exprdata[[i]] <- expr_list
}

# vector with names of all data sets
dataset_vector = str_replace_all(expr_dir_paths, c("/inwosu/Meta_Analysis/Data/race_expression_data/" = "", ".tsv.gz" = ""))

# meta_analysis_object <- list()
# meta_analysis_results <- list()

# Perform meta-analysis
for (i in seq_along(list_all_exprdata)) {
  meta_object <- createObjectMA(listEX = list_all_exprdata[[i]], 
                              listPheno = list_all_metadata[[i]], 
                              namePheno = phenoGroups, 
                              expGroups = phenoCases, 
                              refGroups = phenoControls)
  
  names(meta_object) <- names(list_all_metadata[[i]])
  meta_results <- metaAnalysisDE(meta_object, typeMethod = "REM") 
  
  new_results <- meta_results |>
    filter(FDR < 0.05) |>
    arrange(desc(abs(Com.ES))) |> 
    head(n = 500) |>
    as_tibble(rownames = NA) |> 
    rownames_to_column("Gene")
  
  missing_dataset <- setdiff(dataset_vector, names(meta_object))
  write_tsv(new_results, file.path(meta_results_dir, paste0("meta_results_without_", missing_dataset, ".tsv")))
  
  # meta_analysis_object[[i]] <- meta_object
  # meta_analysis_results[[i]] <- new_results
}

# Draw heatmap
# source("/inwosu/Meta_Analysis/1_run_metaanlysis/functions/draw_Heatmap.R")
# heatmap_values <- draw_Heatmap(objectMA = meta_analysis_object[[i]],
#                                resMA = meta_analysis_results[[i]],
#                                typeMethod = "REM",
#                                scaling = "zscor",
#                                fdrSig = 0.05,
#                                regulation = "all",
#                                numSig = 500,
#                                case = "Black",
#                                control = "White",
#                                title = "Race")

# This section is for meta analysis on the full race data sets.

for (i in seq_along(expr_dir_paths)) {
  metadata <- read_tsv(file = meta_dir_paths[i]) |>
    column_to_rownames(var = "Sample_ID") |>
    filter(race %in% c("Black", "White")) |>
    filter(!is.na(race))
  meta_list[[i]] <- metadata

  expr_data <- read_tsv(file = expr_dir_paths[i]) |>
    dplyr::select(-c("Dataset_ID", "Entrez_Gene_ID", "Chromosome", "Ensembl_Gene_ID", "Gene_Biotype")) |>
    distinct(HGNC_Symbol, .keep_all = TRUE) |>
    inner_join(common_genes) |>
    column_to_rownames(var = "HGNC_Symbol") |>
    na.omit() |>
    as.matrix()
  expr_list[[i]] <- expr_data
}

phenoGroups = rep("race", length(meta_dir_paths))
phenoCases = rep(list("Black"), length(meta_dir_paths))
phenoControls = rep(list("White"), length(meta_dir_paths))

meta_object <- createObjectMA(listEX = expr_list,
                              listPheno = meta_list,
                              namePheno = phenoGroups,
                              expGroups = phenoCases,
                              refGroups = phenoControls)


names(meta_object) <- str_replace_all(expr_dir_paths, c("/inwosu/Meta_Analysis/Data/race_expression_data/" = "", ".tsv.gz" = ""))

resultsMA <- metaAnalysisDE(meta_object, typeMethod = "REM")

new_results <- resultsMA |>
  filter(FDR < 0.05) |>
  arrange(desc(abs(Com.ES))) |> 
  head(n = 50) |>
  as_tibble(rownames = NA) |> 
  rownames_to_column("Gene")

source("/inwosu/Meta_Analysis/2_race_metaanalysis/functions/draw_Heatmap.R")

heatmap_values <- draw_Heatmap(objectMA = meta_object,
                               resMA = new_results,
                               typeMethod = "REM",
                               scaling = "zscor",
                               fdrSig = 0.05,
                               regulation = "all",
                               numSig = 500,
                               case = "Black",
                               control = "White",
                               title = "Race")

 
top_genes = as.data.frame(rownames(heatmap_values)) #check if top_genes is the same as new_result 
# write_tsv(top_genes, file.path("/inwosu/Meta_Analysis/Data/", "race_genes.tsv"))

