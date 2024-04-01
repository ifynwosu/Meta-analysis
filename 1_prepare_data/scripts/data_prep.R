
library(tidyverse)
library(tools)
library(janitor)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 10000)

source("functions/replace_values.R")
source("scripts/variable_standardization.R")
source("scripts/select_files.R")
source("scripts/scale_data.R")

warnings()