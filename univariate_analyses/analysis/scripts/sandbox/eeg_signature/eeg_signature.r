# This script computes both dot product and rho as signatures of EEG-fMRI
# association for each fMRI network on each TR during GradCPT

# This script is quite time / RAM intensive
# See line 120 below for considerations around adjusting the time / RAM
# balance


rm(list=ls())
library(tidyverse)
library(future)
library(future.apply)
library(here)
library(fs)
library(arrow)
library(data.table)


# --- DEFINE FUNCTIONS --- #

compute_tr_metrics <- function(row, cors, networks) {
    # row: All EEG features from one TR from one subject & session during GradCPT
    # cors: Correlations for all EEG features with each fMRI network averaged 
    #       across only the ES runs (across all subjects & sessions)
    
    # returns: A data frame of EEG signatures (dot products and rho) for each fMRI network
    
    # Parse metadata
    row <- unlist(row)
    subject_row <- row['subject']
    session_row <- row['session']
    tr <- row['tr']
    
    # Make a [EEG feature | power value] data frame from the TR
    eeg_features <- sapply(row[!names(row) %in% c('subject', 'session', 'tr')], as.numeric)
    eeg_features <- data.table(feature = names(eeg_features), power = as.numeric(eeg_features))
    
    # Format cors to align with the power data frame computed above
    eeg_cors <- cors[
        subject == subject_row & session == session_row, # Subset
    ][
       , feature := paste(channel, frequency, sep = '_') # Make feature key
          
    ][
        , c('feature', networks), with = FALSE # Select feature key and fmri networks
    ]
    
    # Join into one df to ensure EEG feature alignment
    aligned <- suppressWarnings(eeg_features[eeg_cors, on = .(feature, feature)][, c('feature', 'feature.1') := NULL])
    
    # Compute signature metrics for each fMRI network
    metrics <- data.frame(t(apply(aligned[, colnames(aligned)[colnames(aligned) != 'power'], with = FALSE], 
                                  MARGIN = 2, compute_network_metrics, aligned$power)))
    
    # Format output and return
    colnames(metrics) <- c('dot_product', 'rho')
    metrics <- cbind(data.frame(subject=unname(subject_row), session=unname(session_row), tr=unname(tr), 
                                network=rownames(metrics)), metrics)
    rownames(metrics) <- 1:(nrow(metrics))
    
    return(metrics)
    
}

compute_network_metrics <- function(network, power) {
    # network: a vector of EEG-fMRI correlations across all EEG features 
    #          averaged across ES runs for one subject and session and fMRI network
    # power: a vector of power values across all EEG features for one subject 
    #        and session for one TR of the gradcpt run only
    
    # returns: the dot product and rho between the above vectors
    
    rho <- cor(network, power, method = 'spearman')
    dot_product <- sum(network * power)
    return(c(dot_product, rho))
}

truncate_gradcpt_data <- function(d) {
    # Takes in TR level gradcpt data
    # Keeps only lag 0
    # Updates column names appropriately
    
    stem <- d[,1:3]
    leaf <- d[, 4:ncol(d)]
    
    leaf <- leaf[, grepl('*_0$', colnames(leaf)), with=FALSE]
    new_cols <- unname(sapply(colnames(leaf), FUN = function(x) str_replace(x, '_\\d+$', '')))
    colnames(leaf) <- new_cols
    return(cbind(stem, leaf))
}

# --- SCRIPT STARTS RUNNING HERE --- #

# Parse command line args
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
    setwd(here())
    data_root <- path('analysis/data/original')
} else if (length(args) == 1) {
    data_root <- path(args[1])
} else {
    stop('Usage: Rscript eeg_signature.r path/to/data')
}

print('Importing and formatting data...')

# Import data
cors <- as.data.table(read_feather('analysis/data/correlation_data/correlations_long.feather'))
if (!file.exists(path(data_root, '../correlation_data/correlations_long.feather'))) {
    stop('Must run full preprocessing pipeline and produce correlations_long.feather before running this script.')
}
cors <- as.data.table(read_feather(path(data_root, '../correlation_data/correlations_long.feather')))
if (!file.exists(path(data_root, '../correlation_data/gradcpt_eeg.feather'))) {
    stop('Must run extract_data.py before this script')
}
d <- as.data.table(read_feather(path(data_root, '../correlation_data/gradcpt_eeg.feather')))


# TRUNCATE DATA
# Take only lag 0 from gradcpt TR data
d <- truncate_gradcpt_data(d)
# Get cors for ES runs only by EEG feature (first five lags only)
cors <- cors[run != 'run-001' & lag <= 5, .(cors = mean(cors)), 
             by = .(subject, session, channel, frequency, network)]

networks <- unique(cors$network)

# Spread across networks
cors <- dcast(
    cors,
    subject + session + channel + frequency ~ network,
    value.var = 'cors'
)

# Plan parallel session
# (If you run into RAM issues, try lowering the number of max_workers)
max_workers <- 15
min_workers <- 1
workers <- max(min_workers, min(parallel::detectCores() - 1, max_workers), min_workers)
options(future.globals.maxSize = 16 * 1024^3)

if (interactive()) {
    plan('multisession', workers = 1)
} else {
    plan('multicore', workers = workers)
}

# Split data to ease memory
chunk_size <- 1000
row_chunks <- split(1:nrow(d), ceiling(seq_along(1:nrow(d))/chunk_size))

# Compute signature metrics
print('Computing EEG signature metrics...')
out <- future_lapply(row_chunks, FUN = function(row_idxs) {
    d_sub <- d[row_idxs,]
    cors_sub <- cors[subject %in% unique(d_sub$subject) & session %in% unique(d_sub$session),]
    
    do.call(rbind, lapply(1:nrow(d_sub), function(i) {
        compute_tr_metrics(d_sub[i, ], cors_sub, networks)
    }))
}, future.seed = TRUE)


out <- do.call(rbind, out)

# Save result
out_path <- path(data_root, '../correlation_data/eeg_signatures.csv')
write.csv(out, out_path, row.names=FALSE)

print(paste0('Success! Computed EEG signature data can be found at ', out_path))
