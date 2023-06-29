library(tidyverse)
library(tools)

# Define the file path to the metadata directory
metadata_dir <- "/Data/renamed_metadata/"

# Create the metadata folder if it doesn't exist
if (!dir.exists(metadata_dir)) {
  dir.create(metadata_dir, recursive = TRUE)
}

meta_dir <- "/Data/prelim_metadata"
#meta_dir <- "/inwosu/curated_data/Data/prelim_metadata"
meta_dir_paths <- list.files(meta_dir, full.names = T)
#a <- as.data.frame(meta_dir_paths)

for (i in seq_along(meta_dir_paths)) {
  file = meta_dir_paths[i]
  gseID <- file %>% basename() %>% file_path_sans_ext()
  
  metadata <- read_tsv(file) %>%
    mutate(across(everything(), as.character))
  
  if (gseID == "GSE19615") {
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pr") %>%
      rename(HER2_Status = "her_2") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE19697") {
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pr") %>%
      rename(HER2_Status = "her2") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE20194") {
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pr_status") %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE24185") {
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pr") %>%
      rename(HER2_Status = "her2") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE2990") {
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      replace_ER()
  }
  
  if (gseID == "GSE48390") {
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(HER2_Status = "her2") %>%
      replace_HER2()
  }

  if (gseID == "GSE50948") {
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pr") %>%
      rename(HER2_Status = "her2") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE5460") {
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(HER2_Status = "her2") %>%
      replace_ER() %>%
      replace_HER2()
  }
  
  if (gseID == "GSE58644") {
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(HER2_Status = "her2") %>%
      replace_ER() %>%
      replace_HER2()
  }
  
  if (gseID == "GSE6532_U133A") {
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pgr") %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE6532_U133B") {
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pgr") %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE6532_U133Plus2") {
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pgr") %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE7390") {
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      replace_ER()
  }
  
  if (gseID == "GSE76275") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pr") %>%
      rename(HER2_Status = "her2") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE90521") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pgr") %>%
      rename(HER2_Status = "her2") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE9195") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pgr") %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE11001") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pr") %>%
      rename(HER2_Status = "her2")
  }
  
  if (gseID == "GSE43365") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pr") %>%
      rename(HER2_Status = "her2")%>%
      replace_HR_status()
  }
  
  if (gseID == "GSE61304") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pgr") %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE81538") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_consensus") %>%
      mutate(ER_Status = case_when(ER_Status %in% c("0", "1") ~ "negetive", 
                                   ER_Status == "3" ~ "positive", 
                                   ER_Status == "2" ~ "equivocal",
                                   TRUE ~ ER_Status)) %>%
      rename(PR_Status = "pgr_consensus") %>%
      mutate(PR_Status = case_when(PR_Status %in% c("0", "1") ~ "negetive", 
                                   PR_Status == "3" ~ "positive", 
                                   PR_Status == "2" ~ "equivocal",
                                   TRUE ~ PR_Status)) %>%
      rename(HER2_Status = "her2_consensus") %>%
      mutate(HER2_Status = case_when(HER2_Status == "0" ~ "negetive",
                                     HER2_Status == "1" ~ "positive",
                                     TRUE ~ HER2_Status))
  }
  
  if (gseID == "GSE21653") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_ihc") %>%
      rename(PR_Status = "pr_ihc") %>%
      rename(HER2_Status = "erbb2")%>%
      replace_HR_status()
  }
  
  if (gseID == "METABRIC") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "ER_IHC") %>%
      replace_ER()
  }
  
  if (gseID == "GSE23720") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_ihc") %>%
      rename(PR_Status = "pr_ihc")
  }
  
  if (gseID == "GSE31448") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_ihc") %>%
      rename(PR_Status = "pr_ihc")  %>%
      rename(HER2_Status = "erbb2_ihc_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE17907") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_ihc_status") %>%
      rename(PR_Status = "pr_ihc_status")  %>%
      rename(HER2_Status = "erbb2_ihc_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE17907") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_ihc_status") %>%
      rename(PR_Status = "pr_ihc_status")  %>%
      rename(HER2_Status = "erbb2_ihc_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE18728") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_original") %>%
      rename(PR_Status = "pr")  %>%
      rename(HER2_Status = "her2_summary") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE16391") { 
    new_metadata <- metadata %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HER2() %>%
      mutate(ER_Status = case_when(
        `ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2` == "1" ~ "Positive",
        `ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2` == "2" ~ "Positive"
      ),
      PR_Status = case_when(
        `ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2` == "1" ~ "Positive",
        `ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2` == "2" ~ "Negative")) %>%
      dplyr::select(-`ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2`) %>%
      dplyr::select(Dataset_ID, Sample_ID, Platform_ID, ER_Status, PR_Status, HER2_Status, everything())
  }
  
  if (gseID == "GSE23988") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_positive_vs_negative") %>%
      replace_ER()
  }
  
  if (gseID == "GSE22093") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_positive_vs_negative_by_immunohistochemistry") %>%
      replace_ER()
  }
  
  if (gseID == "GSE18864") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pr_status")  %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE20194") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pr_status")  %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE20271") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pr_status")  %>%
      rename(HER2_Status = "her_2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE20711") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(HER2_Status = "her2_status") %>%
      replace_ER() %>%
      replace_HER2()
  }
  
  if (gseID == "GSE21947") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE31192") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE45255") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pgr_status")  %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE4922_U133A") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE4922_U133B") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE5327") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE5847") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "ER_status") %>%
      rename(HER2_Status = "Her2Neu") %>%
      replace_ER() %>%
      replace_HER2()
  }
  
  if (gseID == "GSE58984") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pgr_status") %>%
      mutate(HER2_Status = "positive") %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE7378") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE96058_HiSeq") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pgr_status")  %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE96058_NextSeq") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pgr_status")  %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE42568") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE47109") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE96058_NextSeq") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pgr_status")  %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE8193") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pr_status")  %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE32518") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status_by_ihc") %>%
      rename(PR_Status = "pr_status_by_ihc")  %>%
      rename(HER2_Status = "her2_status_ihc_and_fish") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE62944_Normal") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status_by_ihc") %>%
      rename(PR_Status = "pr_status_by_ihc")  %>%
      rename(HER2_Status = "her2_status_by_ihc") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE62944_Tumor") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status_by_ihc") %>%
      rename(PR_Status = "pr_status_by_ihc")  %>%
      rename(HER2_Status = "her2_status_by_ihc") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE25055") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status_ihc") %>%
      rename(PR_Status = "pr_status_ihc")  %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE25065") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status_ihc") %>%
      rename(PR_Status = "pr_status_ihc")  %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE32646") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "er_status_ihc") %>%
      rename(PR_Status = "pr_status_ihc")  %>%
      rename(HER2_Status = "her2_status_fish")
  }
  
  if (gseID == "GSE16446") {
    new_metadata <- metadata %>%
      rename(ER_Status = "esr1bimod") %>%
      rename(HER2_Status = "her2fishbin") %>%
      replace_ER() %>%
      replace_HER2()
  }
  
  if (gseID == "GSE17705") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "estrogen_receptor_er_status") %>%
      replace_ER()
  }
  
  if (gseID == "E_TABM_158") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "estrogen_receptor_status") %>%
      rename(PR_Status = "progesterone_receptor_status")  %>%
      rename(HER2_Status = "erb_b2_positive_ihc") %>%
      replace_HER2()
  }
  
  if (gseID == "GSE2603") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "path_er_status") %>%
      rename(PR_Status = "path_pr_status")  %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE29431") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "re") %>%
      rename(PR_Status = "rp")  %>%
      rename(HER2_Status = "her2_ihc") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE10810") { 
    new_metadata <- metadata %>%
      rename(ER_Status = "receptor") %>%
      mutate(ER_Status = case_when(ER_Status == "0" ~ "NA",TRUE ~ ER_Status)) %>%
      replace_ER()
  }
  
  
  write_tsv(new_metadata, file.path(
    metadata_dir,
    paste0(gseID, ".tsv")))
  
  print(paste0("Saved ", gseID, ".tsv to ", metadata_dir))
    
}

#new_metadata %>% count(ER_Status)

