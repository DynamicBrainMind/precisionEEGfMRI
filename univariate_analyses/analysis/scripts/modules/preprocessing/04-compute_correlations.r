# CALCULATE ALL RUN-LEVEL CORRELATIONS
# ------------------------------------ 
# Usage: Rscript 04-compute_correlations path/to/data --bytask
# Will compute correlations by task if you pass the --bytask flag
# This script is very RAM intensive!



# Import libraries
rm(list=ls())
library(tidyverse)
library(data.table)
library(glue)
library(future.apply)
library(here)
library(fs)
library(arrow)
script_root <- path(here(), 'analysis/scripts/modules/preprocessing')
source(path(script_root, '04a-compute_correlations_functions.r'))



# --- START INPUT PARSING / INTERACTIVE SCRIPT --- #
# ------------------------------------------------ #

# --- INPUT PARSING --- # 

args <- commandArgs(trailingOnly=TRUE)

by_task <- any(args %in% c('-t', '--bytask'))
data_root <- args[!args %in% c('-t', '--bytask')]

if (length(data_root) == 0) {
    data_root <- path('analysis/data/original')
    setwd(here())
} else{
    data_root <- path(data_root)
}

# -- ! CONFIGURE HERE ! -- #

# if you are running in rstudio and want a custom data path, 
# set your path to the data by uncommenting below
# data_root <- path('path/to/data')

# --- INTERACTIVE --- # 

# --- IMPORT DATA --- #
if (!file.exists(path(data_root, '../correlation_data/merged_data.feather'))) {
    stop('Need to run 03-make_flat_data.py before this script.')
}
d <- read_feather(path(data_root, '../correlation_data/merged_data.feather'))

# Get channel names
ch_names <- suppressWarnings(readLines(path(data_root, '../correlation_data/ch_names.txt')))


d <- correlation_pipeline(d, script_root, ch_names, by_task)

cache_dir <- path(script_root, 'cache')
# --- WRITE TO FILE --- #
file <- ifelse(by_task, 'correlations_long_bytask.feather', 'correlations_long.feather')
write_feather(d, path(data_root, glue('../correlation_data/{file}')))

print(paste0('Success! Output has been saved to ', path(data_root, glue('../correlation_data/{file}'))))

# Clear cache
dir_delete(cache_dir)
rm(list=ls())
gc()



