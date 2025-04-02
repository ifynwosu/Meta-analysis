
# create vectors of all posible combinations of values
pos_er_values <- c("pos", "1", "P", "ER+", "Positive", "ER pos", "ERpos", "POS", "+", "Positve", "estrogen receptor (ER)+",
                   "+;3", "+;4", "+;5", "+;5-6", "+;6", "+;6-7", "+;7", "+;8")
neg_er_values <- c("neg", "0", "N", "ER-", "Negative", "ER neg", "ERneg", "NEG", "-", "ER negative")

pos_pr_values <- c("pos", "1", "P", "PR+", "Positive", "PgR+", "progesterone receptor (PR)+",
                   "+;3", "+;4", "+;5", "+;5-6", "+;6", "+;6-7", "+;7", "+;8")
neg_pr_values <- c("neg", "0", "N", "PR-", "Negative", "PgR-", "progesterone receptor (PR)-", "-", "-;2")

pos_her2_values <- c("pos", "1", "P", "HER2+", "Positive", "over-expression", "He+", "POS", "pos (3+)", "yes", "3+", "3")
neg_her2_values <- c("neg", "0", "N", "HER2-", "Negative", "normal", "He-", "NEG", "no", "0/1", "0-1", "1+", "1+ (only 1 core pos.)")

# different functions to replace individual hormone receptor status like ER, PR, HER2 and race
replace_race <- function(data) {
    data |>
        mutate(race = case_when(
                                race %in% c("Caucasian", "W", "w", "European American", "WHITE", "white") ~ "White",
                                race %in% c("African American", "B", "b", "BLACK OR AFRICAN AMERICAN", "black") ~ "Black",
                                race %in% c("asian", "ASIAN", "a", "A") ~ "Asian",
                                race %in% c("hispanic", "Mexico Hispanic", "H") ~ "Hispanic",
                                TRUE ~ race))
}

replace_ER <- function(data) {
    data |>
        mutate(ER_status = ifelse(ER_status %in% pos_er_values, "positive", ER_status)) |>
        mutate(ER_status = ifelse(ER_status %in% neg_er_values, "negative", ER_status))
}

replace_PR <- function(data) {
    data |>
        mutate(PR_status = ifelse(PR_status %in% pos_pr_values, "positive", PR_status)) |>
        mutate(PR_status = ifelse(PR_status %in% neg_pr_values, "negative", PR_status))
}

replace_HER2 <- function(data) {
    data |>
        mutate(HER2_status = ifelse(HER2_status %in% pos_her2_values, "positive", HER2_status)) |>
        mutate(HER2_status = ifelse(HER2_status %in% neg_her2_values, "negative", HER2_status))
}

replace_race <- function(data) {
    data |>
        mutate(race = case_when(
                                race %in% c("Caucasian", "W", "w", "European American", "WHITE", "white") ~ "White",
                                race %in% c("African American", "B", "b", "BLACK OR AFRICAN AMERICAN", "black") ~ "Black",
                                race %in% c("asian", "ASIAN", "a", "A") ~ "Asian",
                                race %in% c("hispanic", "Mexico Hispanic", "H") ~ "Hispanic",
                                TRUE ~ race))
}

# use this function if dataset has all 3 hormone recetor status
replace_HR_status <- function(data) {
    data %>%
        mutate(ER_status = case_when(
                                     ER_status %in% pos_er_values ~ "positive",
                                     ER_status %in% neg_er_values ~ "negative",
                                     TRUE ~ ER_status),
        PR_status = case_when(
                              PR_status %in% pos_pr_values ~ "positive",
                              PR_status %in% neg_pr_values ~ "negative",
                              TRUE ~ PR_status),
        HER2_status = case_when(
                                HER2_status %in% pos_her2_values ~ "positive",
                                HER2_status %in% neg_her2_values ~ "negative",
                                TRUE ~ HER2_status))
}
