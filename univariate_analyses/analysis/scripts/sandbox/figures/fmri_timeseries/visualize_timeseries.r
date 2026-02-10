rm(list=ls())
library(tidyverse)
library(ggpubr)
library(here)
library(fs)
library(arrow)
library(scales)
library(RColorBrewer)
setwd(path(here(), 'analysis'))
root <- path('scripts/sandbox/timeseries')

# Import
d <- read_feather('data/merged_data.feather')
d <- d %>% 
    select(subject:dmn_b) %>% 
    filter(subject == 'sub-002', session == 'ses-001', run == 'run-001')


d %>% 
    filter(tr < 50) %>% 
    gather(region, bold, dmn:dmn_b) %>% 
    mutate(contrast = ifelse(region %in% c('dmn', 'dan'), 'DMN-DAN', 'DMNa-DMNb'),
           region = recode(region, `dan` = 'DAN', `dmn` = 'DMN', `dmn_a` = 'DMNa', `dmn_b` = 'DMNb'),
           time = tr * 2) %>% 
    ggplot(aes(x = time, y = bold)) + 
    geom_line(aes(color = region)) + 
    geom_point(size = 1.5, aes(color = region)) + 
    facet_wrap(~contrast, nrow = 2) + 
    scale_color_manual(values = brewer.pal(4, 'Dark2')) + 
    labs(
        x = 'Time (s)',
        y = 'BOLD signal',
        color = ''
    ) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          strip.background = element_rect(fill = NA),
          legend.position = c(.9, .4),
          text = element_text(size = 16),
          legend.title = element_blank())
    

ggsave(path(root, 'timeseries.png'), width = 1920, height = 1080, units = 'px', dpi = 120)



# Show EEG spectrum with fMRI

eeg <- read.csv(path(root, 'eeg_spectrum.csv'))
eeg <- eeg %>% 
    mutate(sample = 1:(nrow(.))) %>% 
    mutate(time = sample / 250) %>% 
    gather(freq, power, starts_with('X')) %>% 
    mutate(freq = as.integer(str_replace(freq, 'X', ''))) 
    
fmri <- read_feather('data/merged_data.feather')
fmri <- fmri %>% 
    select(subject:dmn_b) %>% 
    filter(subject == 'sub-006', session == 'ses-001', run == 'run-002', tr <= 15) %>% 
    mutate(time = tr * 2) %>% 
    select(time, dmn, dan) %>% 
    gather(region, bold, dmn, dan)

d <- left_join(eeg, fmri)

p1 <- eeg %>% 
    ggplot(aes(x = time, y = freq)) + 
    geom_raster(aes(fill = power)) + 
    scale_fill_gradientn(colors = rev(brewer.pal(11, 'RdBu'))) + 
    labs(
        x = 'Time (s)',
        y = 'Frequency (Hz)'
    ) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = 'none',
          text = element_text(size = 16))

rd <- brewer.pal(11, 'RdBu')[3]
bu <- brewer.pal(11, 'RdBu')[9]

p2 <- fmri %>% 
    mutate(region = toupper(region)) %>% 
    ggplot(aes(x = time, y = bold)) + 
    geom_point(size = 1.5, aes(color = region)) + 
    geom_line(aes(color = region)) + 
    labs(
        x = 'Time (s)',
        y = 'BOLD',
        color = ''
    ) + 
    xlim(0, 30) + 
    ylim(-45, 60) + 
    scale_color_manual(values = c('DAN' = bu, 'DMN' = rd)) + 
    theme_bw() +
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = c(.05, .6),
          text = element_text(size = 16))

g <- ggarrange(p2, p1, nrow = 2)

# Screens
ggsave(plot = g, filename = path(root, 'eeg_fmri_timeseries.png'),
       height = 1080, width = 1920, units = 'px', dpi = 120)


