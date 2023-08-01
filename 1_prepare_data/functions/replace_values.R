#make the vectors into variables

pos_er_values = c("pos", "1", "P", "ER+", "Positive", "ER pos", "ERpos", "POS", "+", "Positve", "estrogen receptor (ER)+",
                  "+;3", "+;4", "+;5", "+;5-6", "+;6", "+;6-7", "+;7", "+;8")
neg_er_values = c("neg", "0", "N", "ER-", "Negative", "ER neg", "ERneg", "NEG", "-", "ER negative")

pos_pr_values = c("pos", "1", "P", "PR+", "Positive", "PgR+", "progesterone receptor (PR)+",
                  "+;3", "+;4", "+;5", "+;5-6", "+;6", "+;6-7", "+;7", "+;8")
neg_pr_values = c("neg", "0", "N", "PR-", "Negative", "PgR-", "progesterone receptor (PR)-", "-", "-;2")

pos_her2_values = c("pos", "1", "P", "HER2+", "Positive", "over-expression", "He+", "POS", "pos (3+)", "yes", "3+", "3")
neg_her2_values = c("neg", "0", "N", "HER2-", "Negative", "normal", "He-", "NEG", "no", "0/1", "0-1", "1+", "1+ (only 1 core pos.)")

replace_ER <- function(data) {
  data %>% 
    mutate(ER_Status = ifelse(ER_Status %in% pos_er_values, "positive", ER_Status)) %>% 
    mutate(ER_Status = ifelse(ER_Status %in% neg_er_values, "negative", ER_Status))
}

replace_PR <- function(data) {
  data %>% 
    mutate(PR_Status = ifelse(PR_Status %in% pos_pr_values, "positive", PR_Status)) %>% 
    mutate(PR_Status = ifelse(PR_Status %in% neg_pr_values, "negative", PR_Status))
}

replace_HER2 <- function(data) {
  data %>% 
    mutate(HER2_Status = ifelse(HER2_Status %in% pos_her2_values, "positive", HER2_Status)) %>% 
    mutate(HER2_Status = ifelse(HER2_Status %in% neg_her2_values, "negative", HER2_Status))
}

replace_race <- function(data) {
  data %>%
    mutate(race = case_when(
      race %in% c("Caucasian", "W", "w", "European American", "WHITE", "white") ~ "White", 
      race %in% c("African American", "B", "b", "BLACK OR AFRICAN AMERICAN", "black") ~ "Black",
      TRUE ~ race))
}

replace_HR_status <- function(data) {
  data %>%
    mutate(ER_Status = case_when(
      ER_Status %in% pos_er_values ~ "positive",
      ER_Status %in% neg_er_values ~ "negative",
      TRUE ~ ER_Status),
      PR_Status = case_when(
      PR_Status %in% pos_pr_values ~ "positive",
      PR_Status %in% neg_pr_values ~ "negative",
      TRUE ~ PR_Status),
      HER2_Status = case_when(
      HER2_Status %in% pos_her2_values ~ "positive",
      HER2_Status %in% neg_her2_values ~ "negative",
      TRUE ~ HER2_Status))
}



# replace_value_r <- function(data) {
#   data %>%
#     mutate(ER_Status = case_when(
#       ER_Status %in% c("pos", "1", "P", "ER+") ~ "Positive",
#       ER_Status %in% c("neg", "0", "N", "ER-") ~ "Negative",
#       TRUE ~ ER_Status),
#       PR_Status = case_when(
#         PR_Status %in% c("pos", "1", "P", "PR+") ~ "Positive",
#         PR_Status %in% c("neg", "0", "N", "PR-") ~ "Negative",
#         TRUE ~ PR_Status),
#       HER2_Status = case_when(
#         HER2_Status %in% c("pos", "1", "P", "HER2+") ~ "Positive",
#         HER2_Status %in% c("neg", "0", "N", "HER2-") ~ "Negative",
#         TRUE ~ HER2_Status),
#       race = case_when(
#         race %in% c("Caucasian", "W", "w") ~ "White", 
#         race %in% c("African American", "B", "b") ~ "Black",
#         TRUE ~ race))
# }

# replace_ER <- function(value) {
#   case_when(
#     value %in% c("pos", "1", "P", "ER+") ~ "Positive",
#     value %in% c("neg", "0", "N", "ER-") ~ "Negative",
#     TRUE ~ value
#   )
# }
# 
# replace_PR <- function(value) {
#   case_when(
#     value %in% c("pos", "1", "P", "PR+") ~ "Positive",
#     value %in% c("neg", "0", "N", "PR-") ~ "Negative",
#     TRUE ~ value
#   )
# }
# 
# replace_HER2 <- function(value) {
#   case_when(
#     value %in% c("pos", "1", "P", "HER2+") ~ "Positive",
#     value %in% c("neg", "0", "N", "HER2-") ~ "Negative",
#     TRUE ~ value
#   )
# }
# 
# replace_race <- function(value) {
#   case_when(
#     race %in% c("Caucasian", "W", "w") ~ "White", 
#     race %in% c("African American", "B", "b") ~ "Black",
#     TRUE ~ race
#   )
# }

# new_metadata <- new_metadata %>%
#   mutate(ER_Status = case_when(
#     ER_Status %in% c("pos", "1", "P", "ER+") ~ "Positive",
#     ER_Status %in% c("neg", "0", "N", "ER-") ~ "Negative",
#     TRUE ~ ER_Status
#   ),
#   PR_Status = case_when(
#     PR_Status %in% c("pos", "1", "P", "PR+") ~ "Positive",
#     PR_Status %in% c("neg", "0", "N", "PR-") ~ "Negative",
#     TRUE ~ PR_Status
#   ),
#   HER2_Status = case_when(
#     HER2_Status %in% c("pos (3+)", "1", "P", "HER2+") ~ "Positive",
#     HER2_Status %in% c("neg", "0", "N", "HER2-") ~ "Negative",
#     TRUE ~ HER2_Status
#   ),
# race = case_when(
#   race %in% c("Caucasian", "W", "w") ~ "White",
#   race %in% c("African American", "B", "b") ~ "Black",
#   TRUE ~ race))
