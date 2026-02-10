# --- MAKE ALL CORRELATION VISUALIZATIONS --- #
# --- correlations_long.feather needs to exist prior to running this script
# --- (and possibly correlations_long_bytask.feather)


# --- LIBRARIES --- #
rm(list=ls())
library(arrow)
library(ggridges)
library(tidyverse)
library(eegUtils)
library(glue)
library(data.table)
library(ggpubr)
library(here)
library(fs)
library(psych)
library(scales)
library(RColorBrewer)
library(reticulate)
setwd(path(here()))
root <- path('analysis/scripts/sandbox/figures/correlations')
source(path(root, 'cor_viz_helpers/heat_maps/plot_heat_maps.r')) # plot_heat()
source(path(root, 'cor_viz_helpers/by_task/by_task.r')) # plot_by_task()
source(path(root, 'cor_viz_helpers/topos/plot_topos.r')) # plot_topo()
source(path(root, 'cor_viz_helpers/inter_session/plot_inter_session.r')) # plot_inter_session()
source(path(root, 'cor_viz_helpers/significance/plot_significance.r')) # plot_significance()
data_root <- path(here(), 'analysis/data/original')
text <- 16
size <- 16






# --- !! CUSTOMIZE THESE DIRECTORIES !! --- #

# Where do you want to save your figures?
#fig_save_root <- path(here(), 'analysis/scripts/sandbox/figures/figures_scratch')
fig_save_root <- path('/Users/ak4379/Documents/data/R21_EEG-fMRI/figures')

# Where is your data (if not in analysis/data/original)
data_root <- path('/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/preprocessed')

# ------------------------------------------ #






# --- PRELIMINARY PROCESSING (SHOULDN'T NEED TO EDIT) --- #

if (!dir_exists(fig_save_root)) {
    dir_create(fig_save_root)
}

# Import data
if (file.exists(path(data_root, '../correlation_data/correlations_long.feather'))) {
    d_run <- read_feather(path(data_root, '../correlation_data/correlations_long.feather'))
} else {
    stop('correlations_long.csv is missing. First run 01-format_data.py, then run 03-make_flat_data.py, 
         then run 04-compute_correlations.r')
}

# Get channel coordinates from Python
use_condaenv('eeg-fmri')
py_run_string("
from analysis.scripts.modules.preprocessing.eeg_utils import get_channel_coordinates
import mne
data_root = r.data_root
raw = mne.io.read_raw_eeglab(f'{data_root}/sub-001/ses-001/eeg/sub-001_ses-001_bld001_eeg_Bergen_CWreg_filt_ICA_rej.set')
ch_names = raw.info['ch_names']
ch_pos = get_channel_coordinates(ch_names)
")
ch_pos <- py$ch_pos
ch_names <- py$ch_names

# Adjust the lag var to (s)
d_run <- d_run %>% 
    mutate(lag = lag * 2)

d <- d_run %>% 
    group_by(subject, channel, frequency, lag, network) %>% 
    summarize(cors = mean(cors))



# ========================================= # 


# --- PLOT GENERATION --- #
# (edit as needed)


# --- HEAT MAPS --- #

# Use arg "write_csv='my_file_name'" to save out data (no need for ".csv" extension)
# It will get saved in the fig_save_root defined above

networks <- c('DNa', 'DNb', 'dATNa', 'dATNb')
networks <- c('Yeo7_DN', 'DNa', 'DNb')
p1 <- plot_heat(d, networks, axis_text = 14, by_channels=TRUE, nlags_s=10, nrow = 1,
                write_csv='test')

# saving for spatially-averaged version
ggsave(plot=p1, file=path(fig_save_root, 'heat_maps.png'), 
       width = 2600, height = 850, units = 'px', dpi = 300)

# saving for channel-wise version
ggsave(plot=p1, file=path(fig_save_root, 'heat_maps.png'), 
       width = 1600, height = 3000, units = 'px', dpi = 300)

# --- TOPO PLOTS --- #

#bands <- c('delta', 'theta', 'alpha', 'beta', 'gamma')
#networks <- c('Yeo7_DN', 'DNa', 'DNb')
#networks <- c('LH_DefaultA_PFCm','LH_DefaultA_IPL','RH_DefaultA_Temp')
#networks <- c('LH_DefaultB_PFCd','LH_DefaultB_IPL','RH_DefaultB_Temp')
networks <- c('LH_ContB_PFCmp','LH_ContB_IPL','RH_ContB_Temp')
#bands <- c('theta', 'alpha', 'beta')
bands <- c('alpha')
bands <- c('theta', 'alpha', 'beta1','beta2')
p2 <- plot_topo(d, networks, bands, scales=c(-.04,.11), nlags_s = 10)

ggsave(plot=p2, file=path(fig_save_root, 'topos.png'), 
       width = 1600, height = 2200, units = 'px', dpi = 300)

# Example of how to merge plots as panels in a composite figure
g <- ggarrange(p1, p2, nrow=2)

ggsave(plot=g, file=path(root, 'cor_viz_helpers/heat_maps/dna_through_dan.png'), 
       height = 10, width = 10, units = 'in', dpi = 300)


# --- SIGNIFICANCE PLOTS --- #
# Heat maps with significance to 10 s lags
networks <- unique(d$network) # Plot all networks
bands <- c('delta', 'theta', 'alpha', 'beta', 'gamma')
# Add argument runs=2:4 to visualize only experience sampling data (bytask must be FALSE and d must be d_run)
networks <- c('DNa', 'DNb', 'FPCNa', 'FPCNb', 'dATNa', 'dATNb')
p <- plot_significance(d=d_run, networks=networks, bands=bands, label_middle=FALSE, axis_text = 12,
                       by_channels=FALSE, nlag_s = 10, bytask=FALSE)

# plot by 3-task breakdwon
p <- plot_significance(d=d_task, networks=networks, bands=bands, label_middle=FALSE, axis_text = 12,
                       by_channels=FALSE, nlag_s = 10, bytask=TRUE)

ggsave(plot=p, file=path(fig_save_root, 'significance.png'),
       width = 10, height = 8, units = 'in', dpi = 300)


# --- INTER SESSION PLOT --- #
networks <- c('DNa')
p <- plot_inter_session(d=d_run, networks=networks, bands=bands, label_middle=FALSE, 
                        return_icc = TRUE, axis_text = 10, overall_text=12)

ggsave(plot=p, file=path(fig_save_root, 'inter_session.png'),
       width = 1920, height = 980, units = 'px', dpi = 300)


# --- BY TASK PLOTS --- #

d_task <- read_feather(path(data_root, '../correlation_data/correlations_long_bytask.feather'))
# Adjust lag to seconds
d_task$lag <- d_task$lag * 2
networks <- unique(d_task$network)
p <- plot_by_task(d_task, networks, bands, by_channels=FALSE, nlag_s=10)


ggsave(plot=p, file=path(fig_save_root, 'by_task.png'),
       width = 1920, height = 980, units = 'px', dpi = 150)
































