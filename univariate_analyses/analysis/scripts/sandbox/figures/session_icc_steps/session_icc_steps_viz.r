# Import libraries
rm(list=ls())
library(tidyverse)
library(data.table)
library(psych)
library(paletteer)
library(arrow)
library(fs)
library(here)

# Import data
setwd(here())
data_root <- path('/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/preprocessed')
fig_save_root <- path('/Users/ak4379/Documents/data/R21_EEG-fMRI/figures')

# Add path to data and figure save root
# data_root <- path('path/to/data')
# fig_save_root <- path('path/to/figures)

dpath <- path(data_root, '../correlation_data/iccs.csv')
d <- read.csv(dpath)
#networks <- unique(d$network)
networks <- c('DNa','DNb','FPCNa','FPCNb','dATNa','dATNb')
d$lag <- d$lag * 2 # Convert to seconds
gradient_colors <- paletteer_c('ggthemes::Green', n=100)
#network_colors <- NA
network_colors <- c(rgb(187/255,55/255,56/255),rgb(254/255,147/255,134/255),rgb(79/255,130/255,181/255),rgb(165/255,218/255,244/255),rgb(55/255,119/255,62/255),rgb(206/255,224/255,164/255))



# --- PLOT CONFIG VARS --- #
# --- EDIT THE BELOW --- #
# (you can change any of these before outputting a specific plot and it'll update for that plot)
# (or keep them constant and they'll apply to all corresponding plots)

overall_text <- 18
axis_text <- 16

# Networks to include
networks <- c('DNa', 'DNb','FPCNa','FPCNb','dATNa','dATNb')

# Which color to assign to each network in the overall ICC plot
# Order of colors should match order of 'networks' var
# (hex values work here)
#network_colors <- c('blue', 'green')

# Colors for gradient ICC plot
#gradient_colors <- c('blue', 'green')

# How many lags (in seconds) to include
nlag_s <- 10

# Whether to plot the heatmap by channel (FALSE means by lag)
by_channel <- FALSE


# ----------------------- #


# --- INSPECT TRS BY RUN --- # 
# Full data for TRs in run histograms across subjects
# (This is a big dataset)
d_full <- read_feather(path(data_root, '../correlation_data/merged_data.feather'))

d_full %>% 
    select(subject:tr) %>% 
    group_by(subject, session, run) %>% 
    summarize(tr = max(tr)) %>% 
    mutate(mins = tr * 2 / 60) %>% 
    ggplot(aes(x = mins)) + 
    geom_histogram(fill = 'steelblue', color = 'black') +
    facet_grid(session~run) + 
    labs(
        x = 'Run length (min)',
        y = 'Frequency'
    ) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_rect(fill = NA, color = 'black'),
          text = element_text(size = overall_text),
          axis.text = element_text(size = axis_text))

ggsave(path(fig_save_root, 'mins_per_run.png'), width = 1920, height = 1080, units = 'px', dpi = 150)

rm(d_full)
gc()




# -- VISUALIZE AGGREGATE ICC RESULT -- #
if (!file.exists(path(data_root, '../correlation_data/iccs_full.csv'))) {
    stop('Need to run analysis/scripts/sandbox/session_icc_steps/03-compute_icc_full.r before making this plot.')
}
icc_full <- read.csv(path(data_root, '../correlation_data/iccs_full.csv'))
colnames(icc_full)[colnames(icc_full) == 'network'] <- 'group'

t <- icc_full %>% 
    mutate(lag = lag * 2) %>% # Convert to s
    filter(lag <= nlag_s, group %in% networks) %>% 
    group_by(group) %>% 
    summarize(mean = mean(icc), mx = max(icc), min = min(icc), sd = sd(icc),
              lag = lag[which.max(icc)], 
              channel = channel[which.max(icc)], frequency = frequency[which.max(icc)]) #%>% 

write.csv(t, path(data_root, '../correlation_data/icc_coordinates.csv'), row.names = FALSE)

t <- t %>% 
    select(group, mean, mx) %>% 
    gather(metric, value, mean, mx) %>% 
    mutate(run = ifelse(metric == 'mean', 5, 6), network = 1)
    #mutate(run = ifelse(metric == 'mean', 5, 6), network = 1)

    

p <- d %>% 
    mutate(run = as.integer(str_extract(run_set, '_(\\d)', group = 1))) %>% 
    filter(lag <= nlag_s, network %in% !!networks) %>% 
    group_by(run, network) %>% 
    summarize(icc = mean(icc)) %>% 
    ggplot(aes(x = run, y = icc, group = network)) +
    geom_point(aes(color = network)) + 
    geom_point(data = t, aes(x = run, y = value, color = group)) + 
    geom_line(aes(color = network)) + 
    labs(
        x = 'Number of runs in data',
        y = 'Mean ICC',
        color = 'Network'
    ) + 
    scale_x_continuous(breaks = 1:6, labels = c(1:4, 'Full data\n(mean)', 'Full data\n(max)')) + 
    theme_bw() + 
    theme(axis.ticks = element_blank(),
          panel.grid = element_blank(),
          text = element_text(size = overall_text),
          axis.text = element_text(size = axis_text),
          legend.position = c(.25, .8))

if (!all(is.na(network_colors))) {
    colors <- network_colors
    names(colors) <- networks
    p <- p + scale_color_manual(values = colors)
}

p

ggsave(plot=p, file=path(fig_save_root, 'icc_overall.png'), width = 2500, height = 1080, units = 'px', dpi = 300)




# --- ICC BY FREQ, LAG/CHANNEL, NETWORK, AND NUMBER OF RUNS --- #

y_var <- ifelse(by_channel, 'channel', 'lag')
y_label <- ifelse(by_channel, 'Channel', 'Lag (s)')

if (length(gradient_colors) < 10) {
    gradient_colors <- colorRampPalette(gradient_colors)(100)
}

p <- d %>% 
    mutate(run_label = recode(run_set, `run_set_1` = 'One run',
                              `run_set_2` = 'Two runs',
                              `run_set_3` = 'Three runs',
                              `run_set_4` = 'Four runs'),
           frequency = as.numeric(frequency)) %>% 
    mutate(run_label = factor(run_label, levels = c('One run', 'Two runs', 
                                                    'Three runs', 'Four runs'))) %>% 
    filter(network %in% !!networks, lag <= nlag_s) %>% 
    group_by(frequency, !!sym(y_var), run_label, network) %>% 
    summarize(icc = mean(icc)) %>% 
    ggplot(aes(x = frequency, y = !!sym(y_var))) + 
    geom_tile(aes(fill = icc)) + 
    facet_grid(network~run_label) +
    scale_fill_gradientn(colors = gradient_colors) + 
    labs(
        x = 'Frequency (Hz)',
        y = y_label,
        fill = 'Intraclass correlation\ncoefficient'
    ) + 
    theme_bw() + 
    theme(panel.grid = element_blank(),
          axis.ticks = element_blank(),
          strip.background = element_rect(fill = NA, color = 'black'),
          legend.position = 'bottom',
          legend.text = element_text(angle = 45, hjust = 1),
          text = element_text(size = overall_text),
          axis.text = element_text(size = axis_text))

p

ggsave(plot=p, file=path(fig_save_root, 'icc_freq_lag.png'), width = 1920, height = 1080, units = 'px', dpi = 150)
