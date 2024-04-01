datadir <- "/Data/prelim_metadata"
file_paths <- list.files(datadir, full.names = T)
big_df <- tibble()
total <- 0

for (file in file_paths) {
  df <- read_tsv(file)
  ID <- df$Dataset_ID[1]
  row_num <- nrow(df)
  new_vec <- cbind(ID, row_num)
  big_df <- rbind(big_df, new_vec)
  total <- total + row_num
}

big_df <- rbind(big_df, total)

# GSE1456 = 159
# GSE3494 = 251
# GSE4922 = 289
# GSE6532 = 327