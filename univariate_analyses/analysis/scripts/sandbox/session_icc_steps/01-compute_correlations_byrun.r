# This script computes EEG-fMRI correlations varying how many runs of data are included
# Usage: Rscript 01-compute_correlations_byrun.r optional/path/to/data 3
# Where the integer 3 indicates to include three runs

# This script is configured to compute correlations in parallel only when called from the command line
#   (ie, not in rstudio)

# Due to memory constraints, I suggest running this script four separate times,
# one for each run number


# Import libraries
rm(list=ls())
library(tidyverse)
library(data.table)
library(future)
library(future.apply)
library(glue)
library(psych)
library(here)
library(fs)
library(arrow)

# --- FUNCTION DEFINITIONS --- #
get_correlations <- function(script_root, runs=1:4) {
    # Whether to compute correlations across different sets of runs
    # Saves out to cache
    # (very RAM intensive)
    source(path(here(), 'analysis/scripts/modules/preprocessing/04a-compute_correlations_functions.r'))
    
    print('Importing and formatting data...')
    d <- read_feather(path(data_root, '../correlation_data/merged_data.feather'))
    
    min_trs <- d %>% 
        select(subject:tr) %>% 
        group_by(subject, session, run) %>% 
        summarize(tr = max(tr)) %>% 
        group_by(session, run) %>% 
        summarize(tr_minimum = min(tr)) 
    
    # Filter down to minimum TRs per run
    d <- d %>% 
        inner_join(min_trs) %>% 
        filter(tr <= tr_minimum) %>% 
        select(-tr_minimum)
    
    cache_dir <- path(script_root, 'cache')
    if (!dir.exists(cache_dir)) {
        dir_create(cache_dir)
    }
    
        preprocess_root <- path(here(), 'analysis/scripts/modules/preprocessing')
        ch_names <- suppressWarnings(readLines(path(data_root, '../correlation_data/ch_names.txt')))
        
        lapply(runs, FUN = function(run_int) {
            print(glue('Calculating subset {run_int} of {max(runs)}'))
            print('')
            run_list <- paste0('run-00', 1:run_int)
            d_sub <- d[d$run %in% run_list,]
            cors <- correlation_pipeline(d_sub, preprocess_root, ch_names)
            cors$run_set <- run_int
            write_feather(cors, path(cache_dir, glue('cache_run-00{run_int}.feather')))
            rm(cors)
            gc()
            NULL
        })
}

# -- SCRIPT STARTS EXECUTING HERE -- #

args <- commandArgs(trailingOnly=TRUE)
runs <- NULL

if (length(args) == 0) {
    stop('Need to pass an integer argument indicating how many runs to summarize')
} else {
    data_root <- NULL
    for (arg in args) {
        if (grepl('/', arg)) {
            data_root <- path(arg)
        } else {
            runs <- as.integer(arg)
        }
    }
    if (is.null(data_root)) {
        data_root <- path('analysis/data/original')
        setwd(here())
    }
}

if (is.null(runs)) stop('Need to pass an integer argument indicating how many runs to summarize')

script_root <- path(here(), 'analysis/scripts/sandbox/session_icc_steps')

# Calculate correlations
get_correlations(script_root, runs=runs)

cache_dir <- path(script_root, 'cache')

files <- length(list.files(cache_dir))

if (files < 4) {
    print('Success! You can now run 01-compute_correlations_byrun.r again with the next run set.')
} else {
    print('Success! You can now run 02-compute_icc.r')
}








































