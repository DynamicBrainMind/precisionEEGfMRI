rm(list=ls())
library(arrow)
library(here)
library(data.table)
library(fs)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) > 1) stop('Can only pass one argument (path to preprocessed data.)')

if (length(args) == 0) {
    setwd(here())
    data_root <- path('analysis/data/original')
} else {
    data_root <- args[1]
}

data_dir <- path(data_root, '../correlation_data')
file <- path(data_dir, 'correlations_long_bytask.feather')
if (!file.exists(file)) stop('correlations_long_bytask.feather must exist in correlation_data dir.')

print('Importing data...')
d <- data.table(read_feather(file))

# Write out full csv
#print('Writing out full csv...')
#write.csv(d, path(data_dir, 'correlations_long_bytask.csv'), row.names=FALSE)

# Average collapsing across channel
d <- d[, .(cors = mean(cors)), by = .(subject, session, run, frequency, lag, task, network)]

# Write out csv collapsing across channel
print('Writing out csv collapsed across channels...')
write.csv(d, path(data_dir, 'correlations_long_bytask_across_channel.csv'), row.names=FALSE)

print('All operations completed successfully.')
