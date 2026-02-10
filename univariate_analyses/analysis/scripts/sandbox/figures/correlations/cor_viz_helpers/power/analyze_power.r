# --- POWER ANALYSIS --- #
# Analyze power for both regular and ICC correlations
# Across DNa, DNb, DAN, and FPNa/FPNb

rm(list=ls())
library(tidyverse)
library(here)
library(pwr)
library(fs)
library(glue)
library(pwrFDR)
library(arrow)
setwd(here())
root <- path('analysis/scripts/sandbox/correlations/power')


d <- read_feather(path(root, '../correlations_long.feather'))
d <- d %>% 
    filter(!network %in% c('dATN', 'DMN'), lag <= 10) 

breaks <- c(0, 1, 4, 8, 12, 30, 40)
labels <- c('init', 'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma')
bins <- unique(cut(d$frequency, breaks=breaks))
labels <- paste(labels, bins, sep=' ')

# Standard correlation analysis

cd <- d %>% 
    mutate(band = cut(frequency, breaks, labels)) %>% 
    filter(band != 'init [0, 1)') %>% 
    group_by(subject, session, channel, band, lag, network) %>% 
    summarize(cors = mean(cors)) %>% 
    group_by(subject, band, lag, network) %>% 
    summarize(cors = mean(cors)) %>% 
    group_by(band, lag, network) %>% 
    summarize(cors_ = mean(cors), p  = t.test(cors, mu = 0)$p.value, effect = mean(cors) / sd(cors)) %>% 
    group_by(network) %>% 
    mutate(p_adj = p.adjust(p, method='fdr')) %>% 
    ungroup() 

# Find weakest network
weakest_network <- cd %>% 
    group_by(network) %>% 
    summarize(effect = mean(effect), .groups='drop') %>% 
    filter(abs(effect) == min(abs(effect))) %>% 
    pull(network)

# Find cell with highest effect size in weakest network
weakest_cell <- cd %>% 
    filter(network == !!weakest_network) %>% 
    filter(abs(effect) == max(abs(effect))) 

weakest_effect <- weakest_cell$effect

# Find strongest effect
strongest_cell <- cd[abs(cd$effect) == max(abs(cd$effect)),]    
strongest_effect <- strongest_cell$effect

N_tests <- length(unique(cd$lag)) * length(unique(cd$network))

# Power for low
# (multiple comparisons)
pwrFDR(abs(weakest_effect), r.1 = 27/28, alpha = .05, groups = 1, N.tests=N_tests, average.power = .8)

# (regular)
pwr.t.test(d = weakest_effect, power = .8, type='one.sample', alternative='two.sided')


# -- LOW -- #
# Two-sample t test power calculation 
# 
# n = 4.977854
# d = 2.030593
# sig.level = 0.05
# power = 0.8
# alternative = two.sided
# 
# NOTE: n is number in *each* group

# Power for high
# (multiple comparisons)
pwrFDR(abs(strongest_effect), r.1 = 27/28, alpha = .05, groups = 1, N.tests=N_tests, average.power = .8)

# (regular)
pwr.t.test(d = strongest_effect, power = .8)

# -- HIGH -- #
# Two-sample t test power calculation 
# 
# n = 4.143403
# d = 2.317165
# sig.level = 0.05
# power = 0.8
# alternative = two.sided
# 
# NOTE: n is number in *each* group

# ICC analysis
get_icc <- function(s1, s2) {
    d <- data.frame(s1, s2)
    i <- psych::ICC(d)$results
    r <- i[i$type=='ICC3',]$ICC
    p <- i[i$type=='ICC3',]$p
    return(tibble(r=r, p=p))
}

cd <- d %>% 
    filter(subject != 'sub-023') %>% 
    group_by(subject, session, frequency, lag, network) %>% 
    summarize(cors = mean(cors)) %>% 
    ungroup() %>% 
    mutate(session = str_replace(session, '-', '')) %>% 
    spread(session, cors) %>% 
    group_by(frequency, lag, network) %>% 
    summarize(get_icc(ses001, ses002),
              .groups = 'drop') %>% 
    group_by(network) %>% 
    mutate(p_adj = p.adjust(p, method='fdr')) %>% 
    ungroup() %>% 
    filter(p_adj < .05)

r_range <- range(abs(cd$r))
pwr.r.test(r=r_range[1], power=.8)

# approximate correlation power calculation (arctangh transformation) 
# 
# n = 52.47483
# r = 0.3756776
# sig.level = 0.05
# power = 0.8
# alternative = two.sided

pwr.r.test(r=r_range[2], power=.8)

# approximate correlation power calculation (arctangh transformation) 
# 
# n = 6.626541
# r = 0.8842785
# sig.level = 0.05
# power = 0.8
# alternative = two.sided



















