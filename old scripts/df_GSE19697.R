# Read in data TSV file
metadata <- read_tsv("/inwosu/Meta_Analysis/Data/prelim_metadata/GSE19697.tsv")

metadata <- rename_var(metadata)

# Read in data TSV file
df<- read_tsv("/Data/expression_data/GSE19697.tsv.gz") %>%
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
df <- df %>%  dplyr::select(metadata$SampleID)

# Check if this is in the same order
all.equal(colnames(df), metadata$SampleID)

#add minimum of each row to dataframe to use as a check
df$min <- apply(df, 1, min)

#create vector of minimum row values
df_min <- apply(df, 1, min)

#subtract minimum value from each element in row
df <- sweep(df, 1, df_min,"-")

#add "1" to each value in matrix
df <- sweep(df, 1, 1, "+")

#log transform values
df = log(df)

#filter out NA values as well as other races
metadata_GSE19697 <- metadata %>% 
  filter(Race %in% c("Black", "White")) %>%
  filter(!is.na(Race))

#get sample IDs
common_SampleID <- pull(metadata_GSE19697, SampleID)

#make new dataframe with commom sampleIDs
df_GSE19697 = df[, common_SampleID] %>%
  rownames_to_column(var = "Gene")

rm(df, metadata, common_SampleID, df_min)
