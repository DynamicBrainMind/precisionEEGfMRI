rm(list=ls())
library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(scales)
library(reticulate)
library(fs)
use_condaenv('eeg-fmri')
root <- path('analysis/scripts/sandbox/run_ml')

py_run_string("
import numpy as np
import pickle
import os
import mne
from pyprojroot import here
from pathlib import Path
os.chdir(here())
root = Path('analysis/scripts/sandbox/run_ml')

with open(root / Path('shap_networks.pkl'), 'rb') as file:
    shap_networks = pickle.load(file)

with open(root / Path('truth.pkl'), 'rb') as file:
    truth = pickle.load(file)

with open(root / Path('surrogates.pkl'), 'rb') as file:
    surrogates = pickle.load(file)


fname = Path('sub-001_ses-001_bld001_eeg_Bergen_CWreg_filt_ICA_rej.set')
raw = mne.io.read_raw_eeglab(Path('analysis/data/original/sub-001/ses-001/eeg') / fname)
ch_names = raw.info['ch_names']
")

truth <- py$truth
truth <- sapply(truth, mean)
shaps <- lapply(py$shap_networks, FUN = function(x) colMeans(abs(x)))
surrogates <- py$surrogates
ch_names <- py$ch_names
bands <- c('Delta', 'Theta', 'Alpha', 'Beta', 'Gamma')

# Compile shaps
cols_df <- expand.grid(channel = ch_names, band = bands, lag = 1:6)
cols <- with(cols_df, paste(channel, band, lag, sep='_'))

shap_df <- bind_rows(
    lapply(names(shaps), function(nm) {
        d <- data.frame(col=cols, network=nm, shap = shaps[[nm]])
        d
    })
)


shap_df <- shap_df %>% 
    separate(col, into = c('channel', 'band', 'lag'), sep = '_') %>% 
    mutate(band = factor(band, levels = bands)) %>% 
    group_by(band, lag, network) %>% 
    summarize(shap = mean(shap)) 


text <- 20
# Format surrogates
truth_df <- data.frame(network = names(truth), m = truth)
surrogates_df <- data.frame(do.call(cbind, lapply(surrogates, FUN = function(x) unlist(x))))
p1 <- surrogates_df %>% 
    mutate(iteration = 1:n()) %>% 
    gather(network, surrogate, -iteration) %>% 
    inner_join(truth_df) %>% 
    ggplot(aes(x = surrogate)) + 
    geom_histogram(color = 'black', fill = 'steelblue', bins = 20) +
    geom_point(aes(x = m, y = 50), shape = 8, color = 'steelblue', size = 4) + 
    facet_wrap(~network, nrow=1) + 
    labs(
        x = latex2exp::TeX('$R^2$'),
        y = 'Frequency',
        caption = 'p = .002'
    ) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size = text),
          strip.background = element_rect(fill = 'NA'))

ggsave(filename=path(root, 'surrogate_result.png'), plot = p1, height = 6, width = 10, units = 'in', dpi = 300)


# Plot shaps
min_shap <- min(shap_df$shap)
max_shap <- max(shap_df$shap)

p2 <- shap_df %>% 
    mutate(lag = factor((as.numeric(lag) * 2) - 2)) %>% 
    ggplot(aes(x = band, lag)) + 
    geom_tile(aes(fill = shap)) +
    scale_fill_gradientn(colors = brewer.pal(11, 'Reds')) +   
                         #values = rescale(c(min(shap), 0, max(shap)))) +
    labs(
        x = 'EEG Frequency Band',
        y = 'Lag (s)',
        fill = 'SHAP Feature Importance',
        title = 'Mean absolute value importance'
    ) + 
    facet_wrap(~network, nrow=1) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'bottom',
          text = element_text(size = text),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.text = element_text(angle = 45, hjust = 1, size = 14),
          strip.background = element_rect(fill = NA, color = 'black'))
    
                         
ggsave(path(root, 'SHAP_abs.png'), height = 6, width = 12, units = 'in', dpi = 300)


g <- ggarrange(p1, p2, nrow = 2, labels = c('A.', 'B.'))
g

ggsave(filename = path(root, 'surrogate_shap.png'), plot = g, height = 8, 
       width = 14, units = 'in', dpi = 300)
