
# install.packages("tidyverse", force = TRUE)
# BiocManager::install("RankProd", force = TRUE)

library(tidyverse)
library(RankProd)
library(tools)

# memory.limit(size = 25000)

source("functions/rename_variables.R")
source("scripts/create_dir.R")

meta_dir <- "/Data/prelim_metadata"
meta_dir_paths <- list.files(meta_dir, full.names = T)
# print(meta_dir_paths)
# print("a")
merged_metadata <- NULL

expr_dir  <- "/Data/expression_data"
expr_dir_paths <- list.files(expr_dir, full.names = T)
merged_expr_df <- NULL
expr_data_ncol <- tibble(gseID = "", num_of_cols = 0)
expr_data.origin <- c()

# process metadata
for (i in 1:(length(meta_dir_paths))) {
#for (i in 1:3) {
  num_col <- NULL
  meta_file <- meta_dir_paths[i]
  print(meta_file)
  metadata <- read_tsv(meta_file) 
  metadata <- rename_var(metadata)
  
  # filter out NA values as well as other races
  metadata <- metadata %>%
    filter(race %in% c("Black", "White")) %>%
    filter(!is.na(race))
  
  common_SampleID <- pull(metadata, Sample_ID)
  
  if (is.null(merged_metadata)) {
    merged_metadata <- metadata
  } else {
    merged_metadata <- bind_rows(merged_metadata, metadata)
  }
  
  #process expression file
  expr_file <- expr_dir_paths[i]
  
  expr_data <- read_tsv(expr_file) %>%
    dplyr::select(-c("Dataset_ID","Entrez_Gene_ID", "chromosome_name", "HGNC_Symbol", "Gene_Biotype")) %>%
    distinct(Ensembl_Gene_ID, .keep_all = TRUE) %>%
    na.omit() %>%
    tibble::column_to_rownames("Ensembl_Gene_ID")
  
  #add minimum of each row to dataframe to use as a check
  # expr_data$min_value <- apply(expr_data, 1, min)
  
  # #create vector of minimum row values
  # expr_data_min <- apply(expr_data, 1, min)
  
  # #subtract minimum value from each element in row
  # expr_data <- sweep(expr_data, 1, expr_data_min,"-")
  
  # #add "1" to each value in matrix
  # expr_data <- sweep(expr_data, 1, 1, "+")
  
  # #log transform values
  # expr_data = log(expr_data)
  
  expr_data = expr_data[, common_SampleID] %>%
    rownames_to_column(var = "Ensembl_Gene_ID")
    
  
  if (is.null(merged_expr_df)) {
    merged_expr_df <- expr_data
  } else {
    merged_expr_df <- inner_join(merged_expr_df, expr_data)
  }
  
  num_col <- ncol(expr_data)
  file_name <- meta_file %>% basename() %>% file_path_sans_ext()
  num_col_df <- tibble(gseID = file_name, num_of_cols = num_col)
  expr_data_ncol <- bind_rows(expr_data_ncol, num_col_df)
}

#RankProd requires 2 arguments  
#1: data, a matrix (or data frame) containing the gene expression data. 
#Each row corresponds to a gene, and each column corresponds to a sample.
#2: cl is a vector of length ncol(data) containing the class labels of the samples

#metadata_GSE20271 = metadata_GSE20271[1:10,]
#metadata_GSE20194 = metadata_GSE20194[1:10,]

#expr_data_GSE20271 = expr_data_GSE20271[1:100, 1:10]
#expr_data_GSE20194 = expr_data_GSE20194[1:100, 1:10]

merged_metadata <- merged_metadata %>%
  mutate(race = replace(race, race == "White", 0)) %>%
  mutate(race = replace(race, race == "Black", 1))

merged_expr_df <- merged_expr_df %>% 
  column_to_rownames(var = "Ensembl_Gene_ID")

expr_data_ncol <- expr_data_ncol[-1, ]

for (i in 1:(length(expr_data_ncol$num_of_cols))) {
  a <- expr_data_ncol$num_of_cols[i]
  expr_data.origin <-c(expr_data.origin, rep(i, (a - 1)))
}

expr_data.cl <- as.numeric(pull(merged_metadata, race))
expr_data.gnames <- rownames(merged_expr_df)

start_time <- Sys.time()
RP.adv.out <- RPadvance(merged_expr_df, expr_data.cl, expr_data.origin, logged = TRUE, gene.names = expr_data.gnames)
end_time <- Sys.time()
end_time - start_time

top_genes <- topGene(RP.adv.out, cutoff=0.05, method="pfp", logged=TRUE, logbase=2, gene.names=expr_data.gnames)

#Genes called significant under class1(White) < class2(Black)
up_regulated <- top_genes[["Table1"]]
write.table(up_regulated, file.path(data_dir, "up_regulated.tsv"))

#Genes called significant under class1(White) > class2(Black)
down_regulated <- top_genes[["Table2"]]
write.table(down_regulated, file.path(data_dir, "down_regulated.tsv"))

#length(intersect(rownames(expr_data_GSE20194), rownames(expr_data_GSE19697)))
#browseVignettes("maEndToEnd")
