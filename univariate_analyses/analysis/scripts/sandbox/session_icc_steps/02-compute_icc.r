# This script will compute intraclass correlation coefficients between each EEG feature and
# each fMRI network across 1 - 4 runs included

# Must run 01-compute_correlations_byrun.r before running this script

# Usage: Rscript 02-compute_icc.r optional/path/to/data

# This script is configured to compute ICC in parallel only when called from the command line
#   (ie, not in rstudio)


rm(list=ls())
library(tidyverse)
library(data.table)
library(paletteer)
library(future)
library(future.apply)
library(glue)
library(psych)
library(here)
library(fs)
library(arrow)


get_iccs <- function(script_root) {
    # This function takes in the script root
    # And calculates the ICC across sessions for each set of runs
    # Returns a list with four elements, corresponding to ICCs for each run set
    #   Each list element is a vector of ICCs of length equal to EEG features
    
    if (!dir.exists(path(script_root, 'cache'))) stop('Must run 01-compute_correlations_byrun.r before this script.')
    caches <- list.files(path(script_root, 'cache'), full.names=TRUE)
    print('Importing and formatting data...')
    d <- do.call(rbind, lapply(caches, read_feather))
    d <- data.table(d)
    # Drop subject 3 (no session 1 data)
    d <- d[subject != 'sub-003']
    
    
    d[, network := str_replace(network, '_', '.')]
    d[, feature_set := paste(channel, frequency, lag, network, sep = '_')]
    
    feature_sets <- unique(d$feature_set)
    
    apply_icc <- function(d_run_set) {
        # This function computes ICC in parallel across each EEG feature and network for a 
        # given run set
        
        # Set up parallel
        max_workers <- 1000 # Make this number small if running into memory issues
        workers = max(1, min(parallel::detectCores() - 1, max_workers))
        if (interactive()) {
            plan('multisession', workers = 1)
        } else {
            plan('multicore', workers=workers)
        }
        options(future.globals.maxSize = 16 * 1024^3)
        
        out <- future_sapply(feature_sets, FUN = function(x, d_run_set) {
            # Apply ICC to each EEG feature set and network (in parallel)
            d_sub <- d_run_set[feature_set == x]
            ses_long <- d_sub[, .(cors = mean(cors, na.rm=TRUE)), by = .(subject, session)]
            ses_wide <- dcast(ses_long, subject ~ session, value.var='cors')[, subject := NULL]
            icc <- suppressWarnings(suppressMessages(ICC(ses_wide)$result))
            icc <- icc[icc$type == 'ICC3',]$ICC
            icc
        }, d_run_set = d_run_set)
        
        return(out)
    }
    
    # Apply ICC to each run set
    run_sets <- unique(d$run_set)
    print('Computing ICCs...')
    result <- lapply(run_sets, FUN = function(rs) apply_icc(d[run_set==rs]))
    names(result) <- paste0('run_set_', run_sets)
    
    return(result)
    
}


format_iccs <- function(iccs) {
    # This function takes in the list of ICCs for each run set
    # Returns a concatenated data frame with ICCs labeled by run set and EEG feature
    format_run_set <- function(run_set, label) {
        ds <- data.frame(feature = names(run_set), icc = unname(run_set))
        ds %>% 
            separate(feature, into = c('channel', 'frequency', 'lag', 'network'), sep = '_') %>% 
            mutate(run_set = label)
    }
    
    out <- do.call(rbind, lapply(names(iccs), FUN = function(x) format_run_set(iccs[[x]], label=x)))
    rownames(out) <- NULL
    return(out)
}






# -- SCRIPT STARTS EXECUTING HERE -- #

# Parse input arguments
args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
    setwd(here())
    data_root <- path('analysis/data/original')
} else {
    data_root <- path(args[1])
}

script_root <- path(here(), 'analysis/scripts/sandbox/session_icc_steps')


# Calculate ICCs
iccs <- get_iccs(script_root)
# Format result
d <- format_iccs(iccs)

# Save result
dpath <- path(data_root, '../correlation_data/iccs.csv')
write.csv(d, dpath, row.names=FALSE)
cache_dir <- path(script_root, 'cache')
if (dir_exists(cache_dir)) {
    dir_delete(cache_dir)
}




print('Success! You can now visualize the results with analysis/scripts/sandbox/figures/session_icc_steps/session_icc_steps_viz.r')




