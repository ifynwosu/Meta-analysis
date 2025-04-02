
library(tidyverse)
library(tools)
library(janitor)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10000)

# remove old directories (useful for when scripts don't run to completion)
# old_dir <- c("/Data/temp_metadata/", "/Data/temp_expression_data/", "/Data/analysis_ready_metadata/", "/Data/analysis_ready_expression_data/", "/Data/summary/")

# for (i in (old_dir)) {
#     if (dir.exists(i)) {
#         unlink(i, recursive = TRUE, force = TRUE)
#     }
# }

source("functions/replace_values.R")
source("scripts/variable_standardization.R")
source("scripts/scale_data.R")
source("scripts/select_files.R")
source("scripts/match_data.R")
source("scripts/draw_PCA.R")

warnings()

# total time 60m51.481s