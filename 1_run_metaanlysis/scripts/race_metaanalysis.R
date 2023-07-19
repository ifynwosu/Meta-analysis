
# do we account for batch effects?

# This is where we are going to save the results
meta_results_dir <- "/inwosu/Meta_Analysis/Data/race_meta_results"

if (!dir.exists(meta_results_dir)) {
  dir.create(meta_results_dir, recursive = TRUE)
}

# directory where metadata is stored
meta_dir <- "/inwosu/Meta_Analysis/Data/race_metadata"
meta_dir_paths <- list.files(meta_dir, full.names = T)
meta_list <- list()
#a <- as.data.frame(meta_dir_paths)

# directory where expression data is stored
expr_dir  <- "/inwosu/Meta_Analysis/Data/race_expression_data"
expr_dir_paths <- list.files(expr_dir, full.names = T)
expr_list <- list()

# Define a vector to save combinations for use in generating meta-analysis list
my_vector <- c(length(meta_dir_paths)) 

# Generate all combinations of size (n - 1) because we are leaving one data set out
combinations <- combn(my_vector, (length(meta_dir_paths) - 1))

# create vectors for phenotype data
phenoGroups = rep("race", length(combinations[,1]))
phenoCases = rep(list("Black"), length(combinations[,1]))
phenoControls = rep(list("White"), length(combinations[,1]))

# create different vectors to store values
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
                  column_to_rownames(var = "HGNC_Symbol") |>
                  na.omit() |>
                  as.matrix()
    expr_list[[j]] <- expr_data
  }
  
  list_all_metadata[[i]] <- meta_list
  list_all_exprdata[[i]] <- expr_list
}

meta_analysis_object <- list()
meta_analysis_results <- list()

# create vector to store names of all datasets
dataset_vector = str_replace_all(expr_dir_paths, c("/inwosu/Meta_Analysis/Data/race_expression_data/" = "", ".tsv.gz" = ""))

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

# for (i in seq_along(expr_dir_paths)) {
#   metadata <- read_tsv(file = meta_dir_paths[i]) |>
#     column_to_rownames(var = "Sample_ID") |>
#     filter(race %in% c("Black", "White")) |>
#     filter(!is.na(race))
#   meta_list[[i]] <- metadata
 
#   expr_data <- read_tsv(file = expr_dir_paths[i]) |>
#     dplyr::select(-c("Dataset_ID", "Entrez_Gene_ID", "Chromosome", "Ensembl_Gene_ID", "Gene_Biotype")) |>
#     distinct(HGNC_Symbol, .keep_all = TRUE) |>
#     column_to_rownames(var = "HGNC_Symbol") |>
#     na.omit() |>
#     as.matrix()
#   expr_list[[i]] <- expr_data
# }

# phenoGroups = rep("race", length(meta_dir_paths))
# phenoCases = rep(list("Black"), length(meta_dir_paths))
# phenoControls = rep(list("White"), length(meta_dir_paths))
# 
# meta_object <- createObjectMA(listEX = expr_list, 
#                               listPheno = meta_list, 
#                               namePheno = phenoGroups, 
#                               expGroups = phenoCases, 
#                               refGroups = phenoControls)
# 
# resultsMA <- metaAnalysisDE(meta_object, typeMethod = "REM")
# 
# new_results <- resultsMA |>
#   filter(FDR < 0.05)

# source("/inwosu/Meta_Analysis/1_run_metaanlysis/functions/draw_Heatmap.R")
# heatmap_values <- draw_Heatmap(objectMA = meta_object,
#                                resMA = resultsMA,
#                                typeMethod = "REM",
#                                scaling = "zscor",
#                                fdrSig = 0.05,
#                                regulation = "all",
#                                numSig = 500,
#                                case = "Black",
#                                control = "White",
#                                title = "Race")

# 
# top_genes = as.data.frame(rownames(heatmap_values)) #check if top_genes is the same as new_result 
# write_tsv(top_genes, file.path("/inwosu/Meta_Analysis/Data/", "race_genes.tsv"))
# 
# effects <- calculateES(meta_object)
# head(effects$ES)
