
# Define the file path to the metadata directory
metadata_dir <- "/inwosu/Meta_Analysis/Data/analysis_ready_renamed_metadata/"

# Create the metadata folder if it doesn't exist
if (!dir.exists(metadata_dir)) {
  dir.create(metadata_dir, recursive = TRUE)
}

meta_dir <- "/inwosu/Meta_Analysis/Data/analysis_ready_metadata"
meta_dir_paths <- list.files(meta_dir, full.names = T)

for (i in seq_along(meta_dir_paths)) {
  file = meta_dir_paths[i]
  gseID <- file %>% basename() %>% file_path_sans_ext()
  
  cat("\n")
  print(paste0("Reading in ", gseID, ".tsv"))
  cat("\n")

  metadata <- read_tsv(file) %>%
    mutate(across(everything(), as.character))
  
  if (gseID == "GSE19615") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      mutate(ER_Status = case_when(ER_Status == "pos-low" ~ "equivocal", TRUE ~ ER_Status)) %>%
      rename(PR_Status = "pr") %>%
      mutate(PR_Status = case_when(PR_Status == "pos-low" ~ "equivocal", TRUE ~ PR_Status)) %>%
      rename(HER2_Status = "her_2") %>%
      mutate(HER2_Status = case_when(HER2_Status == "low pos (2+)" ~ "equivocal", TRUE ~ HER2_Status)) %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE19697") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pr") %>%
      rename(HER2_Status = "her2") %>%
      replace_HR_status() %>%
      replace_race()
  }
  
  if (gseID == "GSE20194") {    
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pr_status") %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status() %>%
      replace_race()
  }
  
  if (gseID == "GSE24185") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      mutate(ER_Status = case_when(ER_Status == "U" ~ "NA", TRUE ~ ER_Status)) %>%
      rename(PR_Status = "pr") %>%
      mutate(PR_Status = case_when(PR_Status == "U" ~ "NA", TRUE ~ PR_Status)) %>%
      rename(HER2_Status = "her2") %>%
      mutate(HER2_Status = case_when(HER2_Status == "U" ~ "NA", TRUE ~ HER2_Status)) %>%
      replace_HR_status() %>%
      replace_race()
  }
  
  if (gseID == "GSE2990") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      replace_ER()
  }
  
  if (gseID == "GSE48390") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(HER2_Status = "her2") %>%
      replace_HER2()
  }

  if (gseID == "GSE50948") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pr") %>%
      rename(HER2_Status = "her2") %>%
      replace_HR_status() %>%
      replace_race()
  }
  
  if (gseID == "GSE5460") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(HER2_Status = "her2") %>%
      replace_ER() %>%
      replace_HER2()
  }
  
  if (gseID == "GSE58644") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(HER2_Status = "her2") %>%
      replace_ER() %>%
      replace_HER2()
  }
  
  if (gseID == "GSE6532_U133A") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      mutate(ER_Status = case_when(ER_Status == "?" ~ "equivocal", TRUE ~ ER_Status)) %>%
      rename(PR_Status = "pgr") %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE6532_U133B") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pgr") %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE6532_U133Plus2") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pgr") %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE7390") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      replace_ER()
  }
  
  if (gseID == "GSE76275") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pr") %>%
      rename(HER2_Status = "her2") %>%
      replace_HR_status() %>%
      replace_race()
  }
  
  if (gseID == "GSE90521") {     
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pgr") %>%
      rename(HER2_Status = "her2") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE9195") {     
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pgr") %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE11001") {     
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pr") %>%
      rename(HER2_Status = "her2")
  }
  
  if (gseID == "GSE43365") {    
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      rename(PR_Status = "pr") %>%
      rename(HER2_Status = "her2")%>%
      replace_HR_status()
  }
  
  if (gseID == "GSE61304") {     
    metadata <- metadata %>%
      rename(ER_Status = "er") %>%
      mutate(ER_Status = case_when(ER_Status == "EV" ~ "equivocal", TRUE ~ ER_Status)) %>%
      rename(PR_Status = "pgr") %>%
      mutate(PR_Status = case_when(PR_Status == "EV" ~ "equivocal", TRUE ~ PR_Status)) %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE81538") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_consensus") %>%
      mutate(ER_Status = case_when(ER_Status %in% c("0", "1") ~ "negative", 
                                   ER_Status == "3" ~ "positive", 
                                   ER_Status == "2" ~ "equivocal",
                                   TRUE ~ ER_Status)) %>%
      rename(PR_Status = "pgr_consensus") %>%
      mutate(PR_Status = case_when(PR_Status %in% c("0", "1") ~ "negative", 
                                   PR_Status == "3" ~ "positive", 
                                   PR_Status == "2" ~ "equivocal",
                                   TRUE ~ PR_Status)) %>%
      rename(HER2_Status = "her2_consensus") %>%
      mutate(HER2_Status = case_when(HER2_Status == "0" ~ "negative",
                                     HER2_Status == "1" ~ "positive",
                                     TRUE ~ HER2_Status))
  }
  
  if (gseID == "GSE21653") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_ihc") %>%
      rename(PR_Status = "pr_ihc") %>%
      rename(HER2_Status = "erbb2")%>%
      replace_HR_status()
  }
  
  if (gseID == "METABRIC") {     
    metadata <- metadata %>%
      rename(ER_Status = "ER_IHC") %>%
      replace_ER()
  }
  
  if (gseID == "GSE23720") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_ihc") %>%
      rename(PR_Status = "pr_ihc")
  }
  
  if (gseID == "GSE31448") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_ihc") %>%
      rename(PR_Status = "pr_ihc")  %>%
      rename(HER2_Status = "erbb2_ihc_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE17907") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_ihc_status") %>%
      rename(PR_Status = "pr_ihc_status")  %>%
      rename(HER2_Status = "erbb2_ihc_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE18728") {  
    metadata <- metadata %>%
      rename(ER_Status = "er_original") %>%
      mutate(ER_Status = case_when(ER_Status == "focally pos" ~ "equivocal", TRUE ~ ER_Status)) %>%
      rename(PR_Status = "pr")  %>%
      mutate(PR_Status = case_when(PR_Status == "focally pos" ~ "equivocal", TRUE ~ PR_Status)) %>%
      rename(HER2_Status = "her2_summary") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE16391") {     
    metadata <- metadata %>%
      mutate(ER_Status = case_when(
        `ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2` == "1" ~ "Positive",
        `ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2` == "2" ~ "Positive"
      ),
      PR_Status = case_when(
        `ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2` == "1" ~ "Positive",
        `ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2` == "2" ~ "Negative")) %>%
      dplyr::select(-`ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2`) %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status() %>%
      dplyr::select(Dataset_ID, Sample_ID, Platform_ID, ER_Status, PR_Status, HER2_Status, everything())
  }
  
  if (gseID == "GSE23988") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_positive_vs_negative") %>%
      replace_ER()
  }
  
  if (gseID == "GSE22093") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_positive_vs_negative_by_immunohistochemistry") %>%
      replace_ER()
  }
  
  if (gseID == "GSE18864") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pr_status")  %>%
      rename(HER2_Status = "her2_status") %>%
      mutate(HER2_Status = case_when(HER2_Status == "unk" ~ "NA", TRUE ~ HER2_Status)) %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE20271") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pr_status")  %>%
      rename(HER2_Status = "her_2_status") %>%
      replace_HR_status() %>%
      replace_race()
  }
  
  if (gseID == "GSE20711") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(HER2_Status = "her2_status") %>%
      replace_ER() %>%
      replace_HER2()
  }
  
  if (gseID == "GSE21947") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE31192") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE45255") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pgr_status")  %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE4922_U133A") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      mutate(ER_Status = case_when(ER_Status == "ER?" ~ "equivocal",
                                   TRUE ~ ER_Status)) %>%
      replace_ER()
  }
  
  if (gseID == "GSE4922_U133B") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      mutate(ER_Status = case_when(ER_Status == "ER?" ~ "equivocal",
                                   TRUE ~ ER_Status)) %>%
      replace_ER()
  }
  
  if (gseID == "GSE5327") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE5847") {     
    metadata <- metadata %>%
      rename(ER_Status = "ER_status") %>%
      mutate(ER_Status = case_when(ER_Status == "U" ~ "NA", TRUE ~ ER_Status)) %>%
      rename(HER2_Status = "Her2Neu") %>%
      replace_ER() %>%
      replace_HER2() %>%
      mutate(HER2_Status = case_when(HER2_Status == "U" ~ "NA", TRUE ~ HER2_Status)) %>%
      replace_race()
  }
  
  if (gseID == "GSE58984") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pgr_status") %>%
      mutate(HER2_Status = "positive") %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE7378") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE96058_HiSeq") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pgr_status")  %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE96058_NextSeq") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pgr_status")  %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE42568") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE47109") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      replace_ER()
  }
  
  if (gseID == "GSE8193") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status") %>%
      rename(PR_Status = "pr_status")  %>%
      replace_ER() %>%
      replace_PR()
  }
  
  if (gseID == "GSE32518") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status_by_ihc") %>%
      rename(PR_Status = "pr_status_by_ihc")  %>%
      rename(HER2_Status = "her2_status_ihc_and_fish") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE62944_Normal") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status_by_ihc") %>%
      mutate(ER_Status = case_when(ER_Status == "Indeterminate" ~ "NA", TRUE ~ ER_Status)) %>%
      rename(PR_Status = "pr_status_by_ihc")  %>%
      mutate(PR_Status = case_when(PR_Status == "Indeterminate" ~ "NA", TRUE ~ PR_Status)) %>%
      rename(HER2_Status = "her2_status_by_ihc") %>%
      mutate(across(HER2_Status , ~str_replace(., "Equivocal", "equivocal"))) %>%
      replace_HR_status() %>%
      replace_race()
  }
  
  if (gseID == "GSE62944_Tumor") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status_by_ihc") %>%
      mutate(ER_Status = case_when(ER_Status == "Indeterminate" ~ "NA", TRUE ~ ER_Status)) %>%
      rename(PR_Status = "pr_status_by_ihc")  %>%
      mutate(PR_Status = case_when(PR_Status == "Indeterminate" ~ "NA", TRUE ~ PR_Status)) %>%
      rename(HER2_Status = "her2_status_by_ihc") %>%
      mutate(HER2_Status = case_when(HER2_Status == "Indeterminate" ~ "NA", TRUE ~ HER2_Status)) %>%
      mutate(across(HER2_Status , ~str_replace(., "Equivocal", "equivocal"))) %>%
      replace_HR_status() %>%
      replace_race()
  }
  
  if (gseID == "GSE25055") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status_ihc") %>%
      mutate(ER_Status = case_when(ER_Status == "I" ~ "NA", TRUE ~ ER_Status)) %>%
      rename(PR_Status = "pr_status_ihc")  %>%
      mutate(PR_Status = case_when(PR_Status == "I" ~ "NA", TRUE ~ PR_Status)) %>%
      rename(HER2_Status = "her2_status") %>%
      mutate(HER2_Status = case_when(HER2_Status == "I" ~ "NA", TRUE ~ HER2_Status)) %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE25065") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status_ihc") %>%
      rename(PR_Status = "pr_status_ihc")  %>%
      rename(HER2_Status = "her2_status") %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE32646") {     
    metadata <- metadata %>%
      rename(ER_Status = "er_status_ihc") %>%
      rename(PR_Status = "pr_status_ihc")  %>%
      rename(HER2_Status = "her2_status_fish")
  }
  
  if (gseID == "GSE16446") {    
    metadata <- metadata %>%
      rename(ER_Status = "esr1bimod") %>%
      rename(HER2_Status = "her2fishbin") %>%
      replace_ER() %>%
      replace_HER2()
  }
  
  if (gseID == "GSE17705") {     
    metadata <- metadata %>%
      rename(ER_Status = "estrogen_receptor_er_status") %>%
      replace_ER()
  }
  
  if (gseID == "E_TABM_158") {     
    metadata <- metadata %>%
      rename(ER_Status = "estrogen_receptor_status") %>%
      rename(PR_Status = "progesterone_receptor_status")  %>%
      rename(HER2_Status = "erb_b2_positive_ihc") %>%
      rename(race = "ethnicity") %>%
      replace_HER2() %>%
      replace_race()
  }
  
  if (gseID == "GSE2603") {     
    metadata <- metadata %>%
      rename(ER_Status = "path_er_status") %>%
      rename(PR_Status = "path_pr_status")  %>%
      rename(HER2_Status = "her2_status") %>%
      mutate(HER2_Status = case_when(HER2_Status == "2+" ~ "equivocal",TRUE ~ HER2_Status)) %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE29431") {     
    metadata <- metadata %>%
      rename(ER_Status = "re") %>%
      rename(PR_Status = "rp")  %>%
      rename(HER2_Status = "her2_ihc") %>%
      mutate(HER2_Status = case_when(HER2_Status == "2" ~ "equivocal",TRUE ~ HER2_Status)) %>%
      replace_HR_status()
  }
  
  if (gseID == "GSE10810") {     
    metadata <- metadata %>%
      rename(ER_Status = "receptor") %>%
      mutate(ER_Status = case_when(ER_Status == "0" ~ "NA",TRUE ~ ER_Status)) %>%
      replace_ER()
  }

  if (gseID == "GSE86374") {    
    metadata <- metadata %>%
      rename(race = "population") %>%
      replace_race()
  }
  
  write_tsv(
    metadata, file.path(
    metadata_dir,
    paste0(gseID, ".tsv")))
  
  cat("\n")
  print(paste0("Saved ", gseID, ".tsv to ", metadata_dir))
  cat("\n")

}
