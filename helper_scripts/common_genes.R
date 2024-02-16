# script to test stuff
old_expr_dir  <- "/inwosu/curated_data/Data/analysis_ready_expression_data/"
old_expr_dir_paths <- list.files(old_expr_dir, full.names = T)

old_num_list <- list()

for (i in 1:length(old_expr_dir_paths)) {
  expr_file <- (old_expr_dir_paths[i])
  expr_data <- read_tsv(expr_file)
  my_num <- nrow(expr_data)
  #name_list[[i]] <- expr_data$HGNC_Symbol
  old_num_list[i] <- my_num 
}

expr_dir  <- "/inwosu/Meta_Analysis/Data/results/race/"
expr_dir_paths <- list.files(expr_dir, full.names = T)

num_list <- list()

for (i in 1:length(expr_dir_paths)) {
  expr_file <- (expr_dir_paths[i])
  expr_data <- read_tsv(expr_file)
  num_list[[i]] <- expr_data$Gene
}

common_genes <- tibble(Reduce(intersect, num_list))
names(common_genes) <- "HGNC_Symbol"


