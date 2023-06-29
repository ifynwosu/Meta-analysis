
BiocManager::install(c("BiocManager", "GEOquery",  "SCAN.UPC", "RankProd", "devtools", "tidyverse", "metafor", "rjson"), force = TRUE)


BiocManager::install("biomaRt")

install.packages("tidyverse")
meta_dir <- "/inwosu/Meta_Analysis/test_data/prelim_metadata"
file_paths <- list.files(meta_dir, full.names = T)

expr_data  <- read_tsv("/inwosu/Meta_Analysis/test_data/expression_data/GSE50948.tsv.gz")
expr_data  <- read_tsv("/inwosu/curated_data/Data/analysis_ready_expression_data/METABRIC.tsv.gz")

a  <- read_tsv("/inwosu/Meta_Analysis/Data/Class1(white)_greater_than_Class2(black).tsv")

b  <- read_tsv("/inwosu/Meta_Analysis/Data/Class1(white)_less_Class2(black).tsv")

expr_data <- expr_data %>%
  rename(hgnc_symbol = HGNC)

gene_df <- big_df[order(big_df$Entrez_Gene_ID),]
  
for (file in file_paths) {
  metadata <- read_tsv(file)
  
  file_name <- file %>% basename() %>% file_path_sans_ext()
  # print(file_name)
  # out_file_path_metadata <- paste0(metadata_dir, gseID, ".tsv")
  # out_file_path_raw_metadata <- paste0(raw_metadata_dir, gseID, ".tsv")
}


data = read_tsv("/inwosu/Meta_Analysis/test_data/expression_data/GSE50948.tsv.gz")

collapse_gene_ids = function(x) {
  if (is.na(x))
    return(NA)
  
  return(paste0(sort(unique(x)), collapse=";"))
}

group_by_columns = setdiff(colnames(data), c("Entrez_Gene_ID", "Ensembl_Gene_ID", "HGNC_Symbol"))

collapsed_data = group_by_at(data, all_of(group_by_columns)) %>%
  summarize(Entrez_Gene_ID = collapse_gene_ids(Entrez_Gene_ID),
            Ensembl_Gene_ID = collapse_gene_ids(Ensembl_Gene_ID),
            HGNC_Symbol = collapse_gene_ids(HGNC_Symbol)) %>%
  select(Dataset_ID, Entrez_Gene_ID, Ensembl_Gene_ID, HGNC_Symbol, Gene_Biotype, everything()) %>%
  arrange(Entrez_Gene_ID)



library(tibble)
library(dplyr)

collapse_gene_ids = function(x) {
  x = sort(unique(x))
  return(paste0(x, collapse=","))
}

data = tibble(Entrez = c(50, 50, 50, 51, 52, 53, 53),
              Ensembl = c("ENSG1", "ENSG2", "ENSG3", "ENSG4", "ENSG5", "ENSG6", "ENSG7"),
              HGNC = c("A", "A", "A", "B", "C", "D", "E"),
              Sample1 = c(9.8, 9.8, 9.8, 8.7, 2.2, 3.3, 3.3),
              Sample2 = c(4.4, 4.4, 4.4, 7.1, 0.1, 2.2, 2.2),
              Sample3 = c(6.6, 6.6, 6.6, 9.0, 9.3, 5.5, 5.5))

sample_column_names = setdiff(colnames(data), c("Entrez", "Ensembl", "HGNC"))

data2 = group_by(data, all_of(sample_column_names)) %>%
  summarize(Entrez = collapse_gene_ids(Entrez),
            Ensembl = collapse_gene_ids(Ensembl),
            HGNC = collapse_gene_ids(HGNC)) %>%
  select(Entrez, Ensembl, HGNC, everything())

df <- getFromGEO(gseID)
write_tsv(df, file.path(raw_metadata_dir, paste0(gseID, ".tsv")))

