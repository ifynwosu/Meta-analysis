
rename_var <- function(metadata) {
  
  er_status <- c(`Estrogen Receptor Status` = "er", 
                 `Estrogen Receptor Status` = "er_status")
                 #`Estrogen Receptor Status` = "estrogen_receptor_status")
  pr_status <- c(`Progesterone Receptor Status` = "pr", 
                 `Progesterone Receptor Status` = "pr_status") 
                 #`Progesterone Receptor Status` = "progesterone_receptor_status")
  her2_status <- c(`HER2/Neu Status` = "her2", 
                   `HER2/Neu Status` = "her2_status", 
                   `HER2/Neu Status` = "her_2_status")
  
  
  new_metadata <-  metadata %>%
    rename(any_of(er_status)) %>% 
    rename(any_of(pr_status)) %>% 
    rename(any_of(her2_status)) %>%
    mutate(across(everything(), as.character))
    # mutate(`Estrogen Receptor Status` = as.character(`Estrogen Receptor Status`)) %>%
    # mutate(`Progesterone Receptor Status` = as.character(`Progesterone Receptor Status`)) %>%
    # mutate(`HER2/Neu Status` = as.character(`HER2/Neu Status`))
  
  
    new_metadata <- new_metadata %>% 
        mutate(`Estrogen Receptor Status` = case_match(`Estrogen Receptor Status`,
            c("1", "P", "ER+") ~ "Estrogen Receptor Positive",
            c("0", "N", "ER-") ~ "Estrogen Receptor Negative")) %>%
        mutate(`Progesterone Receptor Status` = case_match(`Progesterone Receptor Status`,
            c("1", "P", "PR+") ~ "Progesterone Receptor Positive",
            c("0", "N", "PR-") ~ "Progesterone Receptor Negative",)) %>%
        mutate(`HER2/Neu Status` = case_match(`HER2/Neu Status`,
            c("1", "P", "HER2+") ~ "HER2/Neu Positive",
            c("0", "N", "HER2-") ~ "HER2/Neu Negative")) %>%
        mutate(race = case_match(race, 
            c("Caucasian", "W", "w") ~ "White",
            c("African American", "B", "b") ~ "Black"))
    
    # cols_to_change <- c("er_status", "pgr_status", "her2_status", "ki67_status")
    
    # metadata_1 <- metadata_1 %>%
    #   mutate_at(all_of(cols_to_change), ~ str_replace(., "0", "negative")) %>%
    #   mutate_at(all_of(cols_to_change), ~ str_replace(., "1", "positive"))

}

