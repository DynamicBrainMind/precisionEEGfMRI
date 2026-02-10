# This script computes an ICC across the whole data (not split by run set)
# For computing the final points needed on the plot showing ICC across runs

rm(list=ls())
library(tidyverse)
library(here)
library(fs)
library(arrow)
library(future)
library(future.apply)
library(psych)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
script_root <- path(here(), 'analysis/scripts/sandbox/session_icc_steps')

if (length(args) == 0) {
    setwd(here())
    data_root <- path('analysis/data/original')
} else if (length(args) == 1) {
    data_root <- path(args[1])
} else {
    stop('Usage: Rscript 03-compute_icc_full.r path/to/data')
}


# Import data
d <- data.table(read_feather(path(data_root, '../correlation_data/correlations_long.feather')))

d <- d[, .(cors = mean(cors)), by = .(subject, session, channel, frequency, lag, network)]

get_icc <- function(v1, v2) {
    t <- data.frame(v1, v2)
    icc <- suppressMessages(ICC(t))
    return(icc$result[3, 2])
}

d <- dcast(d[, session := str_replace(session, '-', '_')], ...~session, value.var='cors')

workers <- parallel::detectCores() - 1
if (interactive()) {
    plan('multisession', workers = 1)
} else {
    plan('multicore', workers=workers)
}

d[, network := str_replace(network, '_', '.')]
d[, condition := paste(channel, frequency, lag, network, sep = '_')]

print('Computing ICCs...')
out <- future_sapply(unique(d$condition), FUN = function(x) {
    t <- d[condition == x, c('ses_001', 'ses_002')]
    out <- get_icc(t$ses_001, t$ses_002)
    out
})


out <- data.frame(condition = names(out), icc = out)
rownames(out) <- NULL

out <- separate(out, condition, into = c('channel', 'frequency', 'lag', 'network'), sep = '_')

write.csv(out, path(data_root, '../correlation_data/iccs_full.csv'), row.names = FALSE)

