rm(list=ls())
library(reticulate)
library(tidyverse)
library(here)
library(RColorBrewer)
library(glue)
library(ggpubr)
library(fs)
fig_save_root <- path(here(), 'analysis/scripts/sandbox/figures/figures_scratch')
overall_text <- 18
script_root <- as.character(path(here(), 'analysis/scripts/sandbox/figures/eeg_timeseries'))

npy_files <- c(path(script_root, 'timeseries.npy'),
               path(script_root, 'power.npy'))

if (!all(file.exists(npy_files))) stop('Must run python extract_timeseries.py path/to/data before this script.')

use_condaenv('eeg-fmri')
py_run_string('
import numpy as np
raw = np.load(f"{r.script_root}/timeseries.npy")
power = np.load(f"{r.script_root}/power.npy")
')

ch_names <- suppressWarnings(readLines(path(script_root, 'ch_names.txt')))

raw <- py$raw
power_in <- py$power

# Convert to long
dimnames(power_in) <- list(
    sample = 1:dim(power_in)[1],
    channel = ch_names,
    freq = 1:40
)

power_in <- as.data.frame.table(power_in, responseName = 'power')


set.seed(42)
random_channel <- sample(ch_names, size = 1)

n_chans <- length(ch_names) %/% 3
random_channels <- sample(ch_names, size = n_chans)
if (!random_channel %in% random_channels) {
    random_channels <- c(random_channels[1:(length(random_channels)-1)], random_channel)
}



# Band plot

# Get frequency bands
power <- power_in
power$freq <- as.numeric(power_in$freq)
breaks <- c(0, 1, 4, 8, 12, 30, 40)
labels <- c('init', 'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma')
labels_simple <- tolower(labels[2:(length(labels))])
bins <- levels(cut(power$freq, breaks=breaks))
labels <- paste(labels, bins, sep=' ')
power <- power %>% 
    mutate(band = as.character(cut(freq, breaks, labels))) %>% 
    mutate(band = ifelse(band == 'init (0,1]', 'Delta (1,4]', band)) %>% 
    mutate(band = factor(band, levels = rev(unique(band)))) %>% 
    filter(channel == random_channel) %>% 
    group_by(sample, band) %>% 
    summarize(power = mean(power)) %>% 
    group_by(band) %>% 
    mutate(power = scale(power)[,1]) %>% 
    ungroup()

p1 <- power %>% 
    mutate(sample = as.numeric(sample)) %>% 
    mutate(time = (sample - min(sample)) / 250) %>% 
    ggplot(aes(x = time, y = power)) + 
    geom_line() + 
    facet_grid(band~.) + 
    labs(
        x = 'Time (s)',
        y = 'Power (normalized)'
    ) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = 18),
          axis.text = element_text(size = 16),
          strip.background = element_rect(fill = NA, color = 'black'))

ggsave(plot=p1, file = path(fig_save_root, 'power_timeseries.png'),
       width = 12, height = 6, units = 'in', dpi = 300)


# TF plot
power <- power_in %>% 
    mutate(sample = as.numeric(sample),
           freq = as.numeric(as.character(freq))) %>% 
    mutate(time = (sample - min(sample)) / 250)

green <- paletteer::paletteer_c('ggthemes::Green', n=100)
rdbu <- brewer.pal(11, 'RdBu')
rdbu <- rev(colorRampPalette(rdbu)(100))

clip_quantile <- .95
time_min <- 0
time_mid <- ceiling(max(power$time) / 2)
time_max <- max(power$time)

p2 <- power %>% 
    filter(channel == random_channel) %>% 
    mutate(power = ifelse(power > quantile(power, clip_quantile), quantile(power, clip_quantile), power)) %>% 
    ggplot(aes(x = time, y = freq)) + 
    geom_raster(aes(fill = power), interpolate = TRUE) +
    labs(
        x = 'Time (s)',
        y = 'Frequency (Hz)',
        fill = 'Power',
        caption = glue('Channel: {random_channel}')
    ) + 
    scale_fill_gradientn(colors = rdbu) + 
    scale_x_continuous(breaks = c(time_min, time_mid, time_max),
                       labels = function(x) as.integer(x)) +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'bottom',
          legend.text = element_text(angle = 45, hjust = 1, size = 10),
          legend.title = element_text(size = 9),
          legend.key.size = unit(.5, 'cm'),
          text = element_text(size = overall_text))
    
    
    
ggsave(plot=p2, file = path(fig_save_root, 'time_freq_spectrum.png'),
       width = 12, height = 6, units = 'in', dpi = 300)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
