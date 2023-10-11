
# Define the file path to the standardized metadata directory
metadata_dir <- "/Data/analysis_ready_renamed_metadata/"

# Create the metadata folder if it doesn't exist
if (!dir.exists(metadata_dir)) {
    dir.create(metadata_dir, recursive = TRUE)
}

old_meta_dir <- "/Data/analysis_ready_metadata"
meta_dir_paths <- list.files(old_meta_dir, full.names = T)

for (i in seq_along(meta_dir_paths)) {
    file = meta_dir_paths[i]
    gseID <- file %>% basename() %>% file_path_sans_ext()
    
    out_file <- paste0(metadata_dir, gseID, ".tsv")
    
    if (file.exists(out_file)) {
        print(paste0(gseID, " has already been processed!"))
    } else {
        cat("\n")
        print(paste0("Reading in ", gseID, ".tsv"))
        cat("\n")

        metadata <- read_tsv(file) %>%
        mutate(across(everything(), as.character))
    
        if (gseID == "GSE19615") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            mutate(ER_status = case_when(ER_status == "pos-low" ~ "equivocal", TRUE ~ ER_status)) %>%
            rename(PR_status = "pr") %>%
            mutate(PR_status = case_when(PR_status == "pos-low" ~ "equivocal", TRUE ~ PR_status)) %>%
            rename(HER2_status = "her_2") %>%
            mutate(HER2_status = case_when(HER2_status == "low pos (2+)" ~ "equivocal", TRUE ~ HER2_status)) %>%
            replace_HR_status()
        }

        if (gseID == "GSE19697") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            rename(PR_status = "pr") %>%
            rename(HER2_status = "her2") %>%
            replace_HR_status() %>%
            replace_race()
        }
        
        if (gseID == "GSE20194") {    
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            rename(PR_status = "pr_status") %>%
            rename(HER2_status = "her2_status") %>%
            replace_HR_status() %>%
            replace_race()
        }

        if (gseID == "GSE24185") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            mutate(ER_status = case_when(ER_status == "U" ~ "NA", TRUE ~ ER_status)) %>%
            rename(PR_status = "pr") %>%
            mutate(PR_status = case_when(PR_status == "U" ~ "NA", TRUE ~ PR_status)) %>%
            rename(HER2_status = "her2") %>%
            mutate(HER2_status = case_when(HER2_status == "U" ~ "NA", TRUE ~ HER2_status)) %>%
            replace_HR_status() %>%
            replace_race()
        }

        if (gseID == "GSE2990") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            replace_ER()
        }

        if (gseID == "GSE48390") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            rename(HER2_status = "her2") %>%
            replace_HER2()
        }
        
        if (gseID == "GSE50948") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            rename(PR_status = "pr") %>%
            rename(HER2_status = "her2") %>%
            replace_HR_status() %>%
            replace_race()
        }
        
        if (gseID == "GSE5460") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            rename(HER2_status = "her2") %>%
            replace_ER() %>%
            replace_HER2()
        }

        if (gseID == "GSE58644") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            rename(HER2_status = "her2") %>%
            replace_ER() %>%
            replace_HER2()
        }

        if (gseID == "GSE6532_U133A") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            mutate(ER_status = case_when(ER_status == "?" ~ "equivocal", TRUE ~ ER_status)) %>%
            rename(PR_status = "pgr") %>%
            replace_ER() %>%
            replace_PR()
        }
        
        if (gseID == "GSE6532_U133B") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            rename(PR_status = "pgr") %>%
            replace_ER() %>%
            replace_PR()
        }

        if (gseID == "GSE6532_U133Plus2") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            rename(PR_status = "pgr") %>%
            replace_ER() %>%
            replace_PR()
        }
        
        if (gseID == "GSE7390") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            replace_ER()
        }
        
        if (gseID == "GSE76275") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            rename(PR_status = "pr") %>%
            rename(HER2_status = "her2") %>%
            replace_HR_status() %>%
            replace_race()
        }
        
        if (gseID == "GSE90521") {     
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            rename(PR_status = "pgr") %>%
            rename(HER2_status = "her2") %>%
            replace_HR_status()
        }
        
        if (gseID == "GSE9195") {     
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            rename(PR_status = "pgr") %>%
            replace_ER() %>%
            replace_PR()
        }
        
        if (gseID == "GSE11001") {     
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            rename(PR_status = "pr") %>%
            rename(HER2_status = "her2")
        }
        
        if (gseID == "GSE43365") {    
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            rename(PR_status = "pr") %>%
            rename(HER2_status = "her2")%>%
            replace_HR_status()
        }
        
        if (gseID == "GSE61304") {     
            metadata <- metadata %>%
            rename(ER_status = "er") %>%
            mutate(ER_status = case_when(ER_status == "EV" ~ "equivocal", TRUE ~ ER_status)) %>%
            rename(PR_status = "pgr") %>%
            mutate(PR_status = case_when(PR_status == "EV" ~ "equivocal", TRUE ~ PR_status)) %>%
            replace_ER() %>%
            replace_PR()
        }
        
        if (gseID == "GSE81538") {     
            metadata <- metadata %>%
            rename(ER_status = "er_consensus") %>%
            mutate(ER_status = case_when(ER_status %in% c("0", "1") ~ "negative", 
                                        ER_status == "3" ~ "positive", 
                                        ER_status == "2" ~ "equivocal",
                                        TRUE ~ ER_status)) %>%
            rename(PR_status = "pgr_consensus") %>%
            mutate(PR_status = case_when(PR_status %in% c("0", "1") ~ "negative", 
                                        PR_status == "3" ~ "positive", 
                                        PR_status == "2" ~ "equivocal",
                                        TRUE ~ PR_status)) %>%
            rename(HER2_status = "her2_consensus") %>%
            mutate(HER2_status = case_when(HER2_status == "0" ~ "negative",
                                            HER2_status == "1" ~ "positive",
                                            TRUE ~ HER2_status))
        }
      
        if (gseID == "GSE21653") {     
            metadata <- metadata %>%
            rename(ER_status = "er_ihc") %>%
            rename(PR_status = "pr_ihc") %>%
            rename(HER2_status = "erbb2")%>%
            replace_HR_status()
        }
        
        if (gseID == "METABRIC") {     
            metadata <- metadata %>%
            rename(ER_status = "ER_IHC") %>%
            replace_ER()
        }
        
        if (gseID == "GSE23720") {     
            metadata <- metadata %>%
            rename(ER_status = "er_ihc") %>%
            rename(PR_status = "pr_ihc")
        }
        
        if (gseID == "GSE31448") {     
            metadata <- metadata %>%
            rename(ER_status = "er_ihc") %>%
            rename(PR_status = "pr_ihc")  %>%
            rename(HER2_status = "erbb2_ihc_status") %>%
            replace_HR_status()
        }
        
        if (gseID == "GSE17907") {     
            metadata <- metadata %>%
            rename(ER_status = "er_ihc_status") %>%
            rename(PR_status = "pr_ihc_status")  %>%
            rename(HER2_status = "erbb2_ihc_status") %>%
            replace_HR_status()
        }
        
        if (gseID == "GSE18728") {  
            metadata <- metadata %>%
            rename(ER_status = "er_original") %>%
            mutate(ER_status = case_when(ER_status == "focally pos" ~ "positive", TRUE ~ ER_status)) %>%
            rename(PR_status = "pr")  %>%
            mutate(PR_status = case_when(PR_status == "focally pos" ~ "positive", TRUE ~ PR_status)) %>%
            rename(HER2_status = "her2_summary") %>%
            replace_HR_status()
        }
      
        if (gseID == "GSE16391") {     
            metadata <- metadata %>%
            mutate(ER_status = case_when(
                `ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2` == "1" ~ "positive",
                `ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2` == "2" ~ "positive"
            ),
            PR_status = case_when(
                `ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2` == "1" ~ "positive",
                `ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2` == "2" ~ "negative")) %>%
            dplyr::select(-`ER_PgR_ER+/PgR+_=_1_ER+/PgR-=_2`) %>%
            rename(HER2_status = "her2_status") %>%
            replace_HR_status() %>%
            dplyr::select(Dataset_ID, Sample_ID, Platform_ID, ER_status, PR_status, HER2_status, everything())
        }
        
        if (gseID == "GSE23988") {     
            metadata <- metadata %>%
            rename(ER_status = "er_positive_vs_negative") %>%
            replace_ER()
        }
        
        if (gseID == "GSE22093") {     
            metadata <- metadata %>%
            rename(ER_status = "er_positive_vs_negative_by_immunohistochemistry") %>%
            replace_ER()
        }
      
        if (gseID == "GSE18864") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            rename(PR_status = "pr_status")  %>%
            rename(HER2_status = "her2_status") %>%
            mutate(HER2_status = case_when(HER2_status == "unk" ~ "NA", TRUE ~ HER2_status)) %>%
            replace_HR_status()
        }
        
        if (gseID == "GSE20271") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            rename(PR_status = "pr_status")  %>%
            rename(HER2_status = "her_2_status") %>%
            replace_HR_status() %>%
            replace_race()
        }
        
        if (gseID == "GSE20711") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            rename(HER2_status = "her2_status") %>%
            replace_ER() %>%
            replace_HER2()
        }
      
        if (gseID == "GSE21947") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            replace_ER()
        }
        
        if (gseID == "GSE31192") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            replace_ER()
        }
        
        if (gseID == "GSE45255") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            rename(PR_status = "pgr_status")  %>%
            rename(HER2_status = "her2_status") %>%
            replace_HR_status()
        }
        
        if (gseID == "GSE4922_U133A") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            mutate(ER_status = case_when(ER_status == "ER?" ~ "equivocal",
                                        TRUE ~ ER_status)) %>%
            replace_ER()
        }
        
        if (gseID == "GSE4922_U133B") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            mutate(ER_status = case_when(ER_status == "ER?" ~ "equivocal",
                                        TRUE ~ ER_status)) %>%
            replace_ER()
        }
        
        if (gseID == "GSE5327") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            replace_ER()
        }
        
        if (gseID == "GSE5847") {     
            metadata <- metadata %>%
            rename(ER_status = "ER_status") %>%
            mutate(ER_status = case_when(ER_status == "U" ~ "NA", TRUE ~ ER_status)) %>%
            rename(HER2_status = "Her2Neu") %>%
            replace_ER() %>%
            replace_HER2() %>%
            mutate(HER2_status = case_when(HER2_status == "U" ~ "NA", TRUE ~ HER2_status)) %>%
            replace_race()
        }
        
        if (gseID == "GSE58984") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            rename(PR_status = "pgr_status") %>%
            mutate(HER2_status = "positive") %>%
            replace_ER() %>%
            replace_PR()
        }
        
        if (gseID == "GSE7378") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            replace_ER()
        }
        
        if (gseID == "GSE96058_HiSeq") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            rename(PR_status = "pgr_status")  %>%
            rename(HER2_status = "her2_status") %>%
            replace_HR_status()
        }
      
        if (gseID == "GSE96058_NextSeq") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            rename(PR_status = "pgr_status")  %>%
            rename(HER2_status = "her2_status") %>%
            replace_HR_status()
        }
        
        if (gseID == "GSE42568") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            replace_ER()
        }
        
        if (gseID == "GSE47109") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            replace_ER()
        }
      
        if (gseID == "GSE8193") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status") %>%
            rename(PR_status = "pr_status")  %>%
            replace_ER() %>%
            replace_PR()
        }
        
        if (gseID == "GSE32518") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status_by_ihc") %>%
            rename(PR_status = "pr_status_by_ihc")  %>%
            rename(HER2_status = "her2_status_ihc_and_fish") %>%
            replace_HR_status()
        }
        
        if (gseID == "GSE62944_Normal") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status_by_ihc") %>%
            mutate(ER_status = case_when(ER_status == "Indeterminate" ~ "NA", TRUE ~ ER_status)) %>%
            rename(PR_status = "pr_status_by_ihc")  %>%
            mutate(PR_status = case_when(PR_status == "Indeterminate" ~ "NA", TRUE ~ PR_status)) %>%
            rename(HER2_status = "her2_status_by_ihc") %>%
            mutate(across(HER2_status , ~str_replace(., "Equivocal", "equivocal"))) %>%
            replace_HR_status() %>%
            replace_race()
        }
      
        if (gseID == "GSE62944_Tumor") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status_by_ihc") %>%
            mutate(ER_status = case_when(ER_status == "Indeterminate" ~ "NA", TRUE ~ ER_status)) %>%
            rename(PR_status = "pr_status_by_ihc")  %>%
            mutate(PR_status = case_when(PR_status == "Indeterminate" ~ "NA", TRUE ~ PR_status)) %>%
            rename(HER2_status = "her2_status_by_ihc") %>%
            mutate(HER2_status = case_when(HER2_status == "Indeterminate" ~ "NA", TRUE ~ HER2_status)) %>%
            mutate(across(HER2_status , ~str_replace(., "Equivocal", "equivocal"))) %>%
            replace_HR_status() %>%
            replace_race()
        }
        
        if (gseID == "GSE25055") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status_ihc") %>%
            mutate(ER_status = case_when(ER_status == "I" ~ "NA", TRUE ~ ER_status)) %>%
            rename(PR_status = "pr_status_ihc")  %>%
            mutate(PR_status = case_when(PR_status == "I" ~ "NA", TRUE ~ PR_status)) %>%
            rename(HER2_status = "her2_status") %>%
            mutate(HER2_status = case_when(HER2_status == "I" ~ "NA", TRUE ~ HER2_status)) %>%
            replace_HR_status()
        }
      
        if (gseID == "GSE25065") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status_ihc") %>%
            rename(PR_status = "pr_status_ihc")  %>%
            rename(HER2_status = "her2_status") %>%
            replace_HR_status()
        }
        
        if (gseID == "GSE32646") {     
            metadata <- metadata %>%
            rename(ER_status = "er_status_ihc") %>%
            rename(PR_status = "pr_status_ihc")  %>%
            rename(HER2_status = "her2_status_fish")
        }
        
        if (gseID == "GSE16446") {    
            metadata <- metadata %>%
            rename(ER_status = "esr1bimod") %>%
            rename(HER2_status = "her2fishbin") %>%
            replace_ER() %>%
            replace_HER2()
        }
        
        if (gseID == "GSE17705") {     
            metadata <- metadata %>%
            rename(ER_status = "estrogen_receptor_er_status") %>%
            replace_ER()
        }
      
        if (gseID == "E_TABM_158") {     
            metadata <- metadata %>%
            rename(ER_status = "estrogen_receptor_status") %>%
            rename(PR_status = "progesterone_receptor_status")  %>%
            rename(HER2_status = "erb_b2_positive_ihc") %>%
            rename(race = "ethnicity") %>%
            replace_HER2() %>%
            replace_race()
        }
        
        if (gseID == "GSE2603") {     
            metadata <- metadata %>%
            rename(ER_status = "path_er_status") %>%
            rename(PR_status = "path_pr_status")  %>%
            rename(HER2_status = "her2_status") %>%
            mutate(HER2_status = case_when(HER2_status == "2+" ~ "equivocal",TRUE ~ HER2_status)) %>%
            replace_HR_status()
        }
        
        if (gseID == "GSE29431") {     
            metadata <- metadata %>%
            rename(ER_status = "re") %>%
            rename(PR_status = "rp")  %>%
            rename(HER2_status = "her2_ihc") %>%
            mutate(HER2_status = case_when(HER2_status == "2" ~ "equivocal",TRUE ~ HER2_status)) %>%
            replace_HR_status()
        }
      
        if (gseID == "GSE10810") {     
            metadata <- metadata %>%
            rename(ER_status = "receptor") %>%
            mutate(ER_status = case_when(ER_status == "0" ~ "NA",TRUE ~ ER_status)) %>%
            replace_ER()
        }
        
        if (gseID == "GSE86374") {    
            metadata <- metadata %>%
            rename(race = "population") %>%
            replace_race()
        }
  
        write_tsv(metadata, file.path(metadata_dir, paste0(gseID, ".tsv")))
        
        cat("\n")
        print(paste0("Saved ", gseID, ".tsv to ", metadata_dir))
        cat("\n")
  }
}
