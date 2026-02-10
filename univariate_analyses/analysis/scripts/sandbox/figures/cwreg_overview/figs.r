rm(list=ls())
library(tidyverse)
library(ggpubr)
library(WaveletComp)
library(zoo)
library(signal)
library(here)
library(fs)
root <- path('analysis/scripts/sandbox/figures_general/cwreg_overview/')
setwd(path(here(), root))

# Import / format
pre_cw <- read.csv('data/pre_cw_timeseries.csv')
pre_cw$dtype <- 'pre_cw'
post_cw <- read.csv('data/post_cw_timeseries.csv')
post_cw$dtype <- 'post_cw'
d <- rbind(pre_cw, post_cw)

d <- d %>% 
    group_by(dtype) %>% 
    mutate(sample = 1:n()) %>% 
    ungroup() %>% 
    gather(channel, voltage, Oz:POz) %>% 
    select(dtype, sample, channel, voltage) 

qu <- read.csv('data/eeg_quality_formatted.csv')


window_size <- 1
fs <- 250
axis_text <- 14
text <- 16

# --- GET CHANNEL WITH HIGHEST MEAN ALPHA --- #

# Visualize channels
# n_chans <- length(unique(d$channel))
# d %>% 
#     dplyr::filter(dtype == 'clean', sample %in% 10000:(10000 + (250))) %>% 
#     ggplot(aes(x = sample, y = voltage)) + 
#     geom_line() + 
#     facet_wrap(~channel, nrow = n_chans)

best_channel <- d %>% 
    dplyr::filter(dtype == 'post_cw') %>% 
    group_by(channel) %>% 
    group_modify(~ {
        wt <- analyze.wavelet(
            my.data = data.frame(signal = .x$voltage),
            my.series = 'signal',
            loess.span = 0,
            dt = 1/250,
            dj = 1/20,
            lowerPeriod = 1/12,
            upperPeriod = 1/8,
            make.pval = FALSE
        )
        
        freqs <- 1 / wt$Period
        alpha_idx <- which(freqs >= 8 & freqs <= 12)
        mean_alpha <- mean(wt$Power[alpha_idx,])
        data.frame(mean_alpha = mean_alpha)
    }) %>% 
    ungroup() %>% 
    print() %>% 
    dplyr::filter(mean_alpha == max(mean_alpha)) %>% 
    pull(channel)

# --- FIND PEAK 10 S ALPHA IN BEST CHANNEL --- #

x <- d[d$dtype=='post_cw' & d$channel==best_channel,]
wt <- analyze.wavelet(
    my.data = x,
    my.series = 'voltage',
    loess.span = 0,
    dt = 1/250,
    dj = 1/20,
    lowerPeriod = 1/40,
    upperPeriod = 1/1,
    make.pval = FALSE
)

tf <- data.frame(t(wt$Power))
colnames(tf) <- 1 / wt$Period
tf <- tf %>% 
    mutate(sample = 1:n()) %>% 
    gather(frequency, power, -sample) %>% 
    mutate(frequency = as.numeric(frequency))


# Bin frequencies
breaks <- c(0, 1, 4, 8, 12, 30, 40)
labels <- c('init', 'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma')
bins <- unique(cut(1:40, breaks=breaks))
labels <- paste(labels, bins, sep=' ')

alpha <- tf %>% 
    mutate(bin = cut(frequency, breaks, labels)) %>% 
    mutate(is_alpha = ifelse(bin == 'Alpha (8,12]', 'alpha', 'non_alpha')) %>% 
    group_by(sample, is_alpha) %>% 
    summarize(power = sum(power)) %>% 
    spread(is_alpha, power) %>% 
    mutate(alpha_adv = alpha / (alpha + non_alpha)) %>% 
    select(sample, alpha_adv) 


ma <- rollmedian(alpha$alpha_adv, window_size * fs, fill = NA)
ma <- data.frame(sample = 1:(length(ma)), alpha = ma)
best_sample <- ma %>% 
    dplyr::filter(sample >= 10000, sample < 150000) %>% 
    #dplyr::filter(sample >= 10000 & sample <= 150000) %>% 
    dplyr::filter(alpha == max(alpha)) %>% 
    print() %>% 
    pull(sample)
    
best_sample <- round(median(best_sample))


start <- best_sample - (window_size/2 * fs)
t_win <- start:(start + window_size*fs)

# --- COMPUTE TIME FREQUENCY --- #

x <- d[d$dtype=='post_cw' & d$channel == best_channel,]$voltage
post_cw_waves <- analyze.wavelet(
    my.data = data.frame(signal=x),
    my.series = 'signal',
    loess.span = 0,
    dt = 1/fs,
    dj = 1/20,
    lowerPeriod = 1/40,
    upperPeriod = 1/1,
    make.pval = FALSE
)
post_cw_mat <- post_cw_waves$Power

x <- d[d$dtype=='pre_cw' & d$channel == best_channel,]$voltage
pre_cw_waves <- analyze.wavelet(
    my.data = data.frame(signal=x),
    my.series = 'signal',
    loess.span = 0,
    dt = 1/fs,
    dj = 1/20,
    lowerPeriod = 1/40,
    upperPeriod = 1/1,
    make.pval = FALSE
)
pre_cw_mat <- pre_cw_waves$Power

# format
post_cw_df <- data.frame(t(post_cw_mat))
colnames(post_cw_df) <- as.character(post_cw_waves$Period)
post_cw_df$dtype <- 'post_cw'
pre_cw_df <- data.frame(t(pre_cw_mat))
colnames(pre_cw_df) <- as.character(pre_cw_waves$Period)
pre_cw_df$dtype <- 'pre_cw'

waves <- rbind(post_cw_df, pre_cw_df)
waves <- waves %>% 
    group_by(dtype) %>% 
    mutate(sample = 1:n()) %>% 
    ungroup() %>% 
    gather(period, power, -sample, -dtype) %>% 
    mutate(period = as.numeric(period)) %>% 
    mutate(frequency = 1 / period) %>% 
    select(-period)


# --- FIGURE: INTERLEAVED TIME DOMAIN - TIME-FREQ DOMAIN ---  #
voltage_range <- range(d[d$dtype=='pre_cw' & d$channel == best_channel & d$sample %in% t_win,]$voltage)
voltage_min <- voltage_range[1]
voltage_max <- voltage_range[2]

p1 <- d %>% 
    dplyr::filter(dtype == 'pre_cw', channel == best_channel, sample %in% t_win) %>% 
    mutate(time = sample / fs) %>% 
    mutate(time = time - min(time)) %>% 
    ggplot(aes(x = time, y = voltage)) + 
    geom_line() + 
    labs(
        x = 'Time (s)',
        y = 'EEG potential (V)',
        title = 'Before CWL regression (after MR gradient correction)'
    ) + 
    ylim(voltage_min, voltage_max) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = text),
          axis.text = element_text(size = axis_text))
        
max_p <- ceiling(max(waves[waves$sample %in% t_win & waves$dtype=='post_cw',]$power))
p2 <- waves %>% 
    dplyr::filter(dtype=='pre_cw', sample %in% t_win) %>% 
    mutate(time = sample / fs) %>% 
    mutate(time = time - min(time)) %>% 
    ggplot(aes(x = time, y = frequency, fill = power)) + 
    geom_raster(interpolate=TRUE) + 
    scale_y_log10() + 
    scale_fill_viridis_c(option = 'plasma', limits = c(0, max_p), breaks = seq(0, max_p, length.out=2)) + 
    labs(
        x = 'Time (s)',
        y = 'Frequency (Hz)',
        fill = 'Power'
    ) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = text),
          axis.text = element_text(size = axis_text),
          legend.position = 'none')

p3 <- d %>% 
    dplyr::filter(dtype == 'post_cw', channel == best_channel, sample %in% t_win) %>% 
    mutate(time = sample / fs) %>% 
    mutate(time = time - min(time)) %>% 
    ggplot(aes(x = time, y = voltage)) + 
    geom_line() + 
    labs(
        x = 'Time (s)',
        y = 'EEG potential (V)',
        title = 'After CWL regression'
    ) + 
    ylim(voltage_min, voltage_max) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = text),
          axis.text = element_text(size = axis_text))

p4 <- waves %>% 
    dplyr::filter(dtype=='post_cw', sample %in% t_win) %>% 
    mutate(time = sample / fs) %>% 
    mutate(time = time - min(time)) %>% 
    ggplot(aes(x = time, y = frequency, fill = power)) + 
    geom_raster(interpolate=TRUE) + 
    scale_y_log10() + 
    scale_fill_viridis_c(option = 'plasma', limits = c(0, max_p), breaks = seq(0, max_p, length.out=2)) + 
    labs(
        x = 'Time (s)',
        y = 'Frequency (Hz)',
        fill = 'Power',
        caption = glue::glue('Max power: {max_p}')
    ) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = text),
          axis.text = element_text(size = axis_text),
          legend.position = 'none')

g1 <- ggarrange(ggarrange(p1, p2, nrow = 2), ggarrange(p3, p4, nrow = 2), 
                nrow = 2, labels = c('A.', 'B.'))

ggsave(filename=glue::glue('r01_figure_window{window_size}.png'), plot = g1, height = 14, width = 10, units = 'in', dpi = 300)

# --- BY FREQ BAND --- #

# Bin frequencies
breaks <- c(0, 1, 4, 8, 12, 30, 40)
labels <- c('init', 'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma')
bins <- unique(cut(1:40, breaks=breaks))
labels <- paste(labels, bins, sep=' ')


p2a <- waves %>% 
    dplyr::filter(dtype=='pre_cw', sample %in% t_win) %>% 
    mutate(time = sample / fs, bin = cut(frequency, breaks, labels)) %>% 
    mutate(time = time - min(time),
           bin = factor(bin, levels = c('Delta (1,4]', 'Theta (4,8]', 'Alpha (8,12]', 'Beta (12,30]', 'Gamma (30,40]'))) %>% 
    group_by(bin, time) %>% 
    summarize(power = mean(power)) %>% 
    ggplot(aes(x = time, y = bin, fill = power)) + 
    geom_raster(interpolate=FALSE) + 
    #scale_y_log10() + 
    scale_fill_viridis_c(option = 'plasma', limits = c(0, max_p), breaks = seq(0, max_p, length.out=2)) + 
    labs(
        x = 'Time (s)',
        y = 'Frequency band (Hz)',
        fill = 'Power'
    ) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = 16),
          axis.text = element_text(size = axis_text),
          legend.position = 'none')

p4a <- waves %>% 
    dplyr::filter(dtype=='post_cw', sample %in% t_win) %>% 
    mutate(time = sample / fs, bin = cut(frequency, breaks, labels)) %>% 
    mutate(time = time - min(time),
           bin = factor(bin, levels = c('Delta (1,4]', 'Theta (4,8]', 'Alpha (8,12]', 'Beta (12,30]', 'Gamma (30,40]'))) %>% 
    group_by(bin, time) %>% 
    summarize(power = mean(power)) %>% 
    ggplot(aes(x = time, y = bin, fill = power)) + 
    geom_raster(interpolate=FALSE) + 
    #scale_y_log10() + 
    scale_fill_viridis_c(option = 'plasma', limits = c(0, max_p), breaks = seq(0, max_p, length.out=2)) + 
    labs(
        x = 'Time (s)',
        y = 'Frequency band (Hz)',
        fill = 'Power'
    ) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = 16),
          axis.text = element_text(size = axis_text),
          legend.position = 'none')


g2 <- ggarrange(ggarrange(p1, p2a, nrow = 2), ggarrange(p3, p4a, nrow = 2), 
                nrow = 2, labels = c('A.', 'B.'))

ggsave(filename=glue::glue('r01_figure_a_window{window_size}.png'), 
       plot = g2, height = 14, width = 10, units = 'in', dpi = 300)
















# --- HISTOGRAMS --- #

p5 <- qu %>% 
    gather(metric, value, good_chans:brain_ica) %>% 
    group_by(subject, session, metric) %>% 
    summarize(value = mean(value)) %>% 
    mutate(metric = recode(metric, `brain_ica` = 'ICA components\nlabeled brain',
                           `good_chans` = 'Good channels',
                           `good_data` = 'Good data')) %>% 
    ggplot(aes(x = value)) + 
    geom_histogram(fill = '#0077BB', color = 'white', bins = 20) + 
    labs(
        x = 'Proportion',
        y = 'Frequency'
    ) + 
    facet_wrap(~metric, nrow = 1, scales = 'free_y')  +
    xlim(0, 1) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          strip.background = element_rect(fill = NA),
          text = element_text(size = 16),
          axis.text = element_text(size = axis_text))
    

g3 <- ggarrange(ggarrange(p1, p2, nrow=2), ggarrange(p3, p4, nrow=2), p5, nrow=3,
                labels = c('A.', 'B.', 'C.'))

ggsave(filename='figure_for_lotus.png', plot = g3, height = 14, width = 10, units = 'in', dpi = 300)


