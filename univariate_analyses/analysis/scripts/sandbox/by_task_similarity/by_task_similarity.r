# Import libraries
rm(list = ls())
library(tidyverse)
library(here)
library(data.table)
library(psych)
library(fs)
library(arrow)


# Parse input args
args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
    setwd(here())
    data_root <- path('analysis/data/original')
} else if (length(args) == 1) {
    data_root <- path(args[1])
} else {
    stop('Ambiguous input arguments.')
}

# Import data
d <- data.table(read_feather(path(data_root, '../correlation_data/correlations_long_bytask.feather')))

# Average across task & EEG feature
d <- d[, .(cors = mean(cors)), by = .(subject, session, channel, frequency, lag, task, network)]

# Fisher transform
d[, cors := fisherz(cors)]

# Yank subject 3
d <- d[subject != 'sub-003']

# Spread across task
d <- dcast(d, ... ~ task, value.var='cors')


# Summarize similarity across tasks for each subject & network
out <- d[, .(gradcpt_rest = cor(gradcpt, rest, method='spearman'),
      gradcpt_tp = cor(gradcpt, tp, method = 'spearman'),
      rest_tp = cor(rest, tp, method = 'spearman')),
  by = .(subject, network)]

t <- melt(out, id.vars = c('subject', 'network'),
          measure.vars = c('gradcpt_rest', 'gradcpt_tp', 'rest_tp'),
          variable.name = 'contrast',
          value.name = 'similarity')

group <- t[, .(m = mean(similarity), sd = sd(similarity)), by = contrast]


write.csv(out, path(data_root, '../correlation_data/task_similarity_subject_network.csv'), row.names=FALSE)
write.csv(group, path(data_root, '../correlation_data/task_similarity_group.csv'), row.names=FALSE)
