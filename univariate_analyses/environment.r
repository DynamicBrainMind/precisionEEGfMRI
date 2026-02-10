install.packages(c('ggridges', 'tidyverse', 'glue', 'data.table', 'ggpubr', 'here',
                   'fs', 'psych', 'scales', 'RColorBrewer', 'reticulate', 'remotes',
                   'latex2exp', 'paletteer'))




# Look into these dependencies if you have trouble installing arrow:
# (note these are for Linux; YMMV on Mac)
    # CURL: install libcurl-devel (rpm) or libcurl4-openssl-dev (deb)
    # OpenSSL >= 1.0.2: install openssl-devel (rpm) or libssl-dev (deb)

source("https://raw.githubusercontent.com/apache/arrow/main/r/R/install-arrow.R")
install_arrow()


# eegUtils
# (this package tends to be very painful; see eeg_utils_troubleshooting.txt)
remotes::install_github('craddm/eegUtils')
