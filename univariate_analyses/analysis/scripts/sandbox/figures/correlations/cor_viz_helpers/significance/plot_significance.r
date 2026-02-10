# -- FREQ X LAG WITH SIGNIFICANCE -- #

plot_significance <- function(d, networks, bands, scales=NA, label_middle=TRUE,
                              overall_text=25, axis_text=18, legend_text=16,
                              y_label=NA, x_label=NA, by_channels=FALSE, title='', 
                              colors=NA, nlag_s=NA, nrow=NA, bytask=FALSE, runs=NA,
                              write_csv=NA) {
    #' Plot Frequency–Lag Heatmap with Significance Markers
    #'
    #' Creates a faceted heatmap of mean correlations (`cors`) across frequency bands and time lags
    #' for specified networks, highlighting statistically significant effects (FDR-corrected p < .05).
    #'
    #' @param d A data frame containing columns `subject`, `frequency`, `lag`, `network`, and `cors`.
    #' @param networks A character vector of network names to include in the plot.
    #' @param bands A character vector of frequency bands to plot (e.g., "delta", "theta", "alpha", "beta", "gamma").
    #' @param scales Optional numeric vector `c(min, max)` to fix the color scale limits. Defaults to automatic scaling.
    #' @param label_middle Logical; if `TRUE`, labels the color scale at zero. Default is `TRUE`.
    #' @param overall_text Numeric, base font size for all text. Default is 25.
    #' @param axis_text Numeric, font size for axis labels. Default is 18.
    #' @param legend.text Numeric, font size for legend text. Default is 16.
    #'
    #' @return A `ggplot` object representing the faceted heatmap, with gold markers indicating significant bins.
    #'
    #' @details
    #' - Frequencies are binned into standard EEG bands: Delta, Theta, Alpha, Beta, Gamma.
    #' - Within each network and lag, one-sample t-tests test whether mean `cors` differs from zero.
    #' - p-values are adjusted for multiple comparisons using the FDR method.
    #' - Gold asterisks indicate significant bins (FDR-adjusted p < .05).
    #'
    #' @examples
    #' plot_significance(df, networks = c("DMN", "DAN"), bands = c("alpha", "beta"))
    
    if (bytask & !'task' %in% colnames(d)) stop('Must pass by task correlation data frame as data argument.')
    if (!all(is.na(runs)) & !'run' %in% colnames(d)) stop('Must pass run-level correlation data frame as data argument.')
    if (!all(is.na(runs)) & bytask) stop('Cannot specify runs while plotting across tasks.')
    
    if (!all(is.na(runs))) runs <- sapply(runs, FUN = function(x) paste0('run-00', x))
    y_var <- ifelse(by_channels, 'channel', 'lag')
    y_label <- ifelse(!is.na(x_label), x_label,
                      ifelse(y_var == 'channel', 'Channel', 'Lag (s)'))
    x_label <- ifelse(is.na(x_label), 'Frequency (Hz)', x_label)
    
    colors <- case_when(
        is.na(colors) ~ rev(brewer.pal(11, 'RdBu')),
        length(colors) == 2 ~ colorRampPalette(c(colors[1], 'white', colors[2]))(11),
        length(colors) > 2 ~ colorRampPalette(colors)(11),
        .default = rev(brewer.pal(11, 'RdBu'))
    )
    
    
	# Get frequency bands
	breaks <- c(0, 1, 4, 8, 12, 30, 40)
	labels <- c('init', 'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma')
	labels_simple <- tolower(labels[2:(length(labels))])
	bins <- levels(cut(d$frequency, breaks=breaks))
	labels <- paste(labels, bins, sep=' ')
	
	# Determine which bands to grab
	bands <- tolower(bands)
	if (!all(bands %in% labels_simple)) {
	    stop(glue('Labels must be one of {paste(labels_simple, collapse = ", ")}'))
	}
	
	bands_grab <- labels[2:(length(labels))][which(labels_simple %in% bands)]
	
	# Generate significance markers
	if (!is.na(nlag_s)) d <- d[d$lag <= nlag_s,]
	if (!bytask) {
	    if (!all(is.na(runs))) d <- d[d$run %in% runs,]
    	ps <- d %>% 
    		mutate(bin = cut(frequency, breaks, labels)) %>% 
    		filter(bin %in% bands_grab, 
    			   network %in% !!networks) %>% 
    		group_by(subject, !!sym(y_var), network, bin) %>% 
    		summarize(cors = mean(cors)) %>% 
    	    ungroup() %>% 
    	    mutate(cors = fisherz(cors)) %>% # Apply Fisher z transform
    		group_by(!!sym(y_var), network, bin) %>% 
    		summarize(p = t.test(cors, mu = 0)$p.value) %>% 
    	    group_by(network) %>% 
    		mutate(p_adj = p.adjust(p, method='fdr'),
    			   network = factor(network, levels = !!networks)) %>% 
    		filter(p_adj < .05) 
    		
    	# Generate averages to show in heat map
    	sd <- d %>% 
    		mutate(bin = cut(frequency, breaks, labels)) %>% 
    		filter(bin %in% bands_grab,
    			   network %in% !!networks) %>% 
    		group_by(subject, bin, !!sym(y_var), network) %>% 
    		summarize(cors = mean(cors)) 
    	pd <- sd %>% 
    		group_by(bin, !!sym(y_var), network) %>% 
    		summarize(cors = mean(cors)) %>% 
    		mutate(network = factor(network, levels = networks),
    		       cors = fisherz(cors))
	} else {
	    
    	ps <- d %>% 
    		mutate(bin = cut(frequency, breaks, labels)) %>% 
    		filter(bin %in% bands_grab,
    			   network %in% !!networks) %>% 
    		group_by(subject, !!sym(y_var), network, bin, task) %>% 
    		summarize(cors = mean(cors)) %>% 
    	    ungroup() %>% 
    	    mutate(cors = fisherz(cors)) %>% # Apply Fisher z transform
    		group_by(!!sym(y_var), network, bin, task) %>% 
    		summarize(p = t.test(cors, mu = 0)$p.value) %>% 
    	    group_by(network, task) %>% 
    		mutate(p_adj = p.adjust(p, method='fdr'),
    			   network = factor(network, levels = !!networks)) %>% 
    		filter(p_adj < .05) 
    		
    	# Generate averages to show in heat map
    	sd <- d %>% 
    		mutate(bin = cut(frequency, breaks, labels)) %>% 
    		filter(bin %in% bands_grab, 
    			   network %in% !!networks) %>% 
    		group_by(subject, bin, !!sym(y_var), network, task) %>% 
    		summarize(cors = mean(cors)) 
    	pd <- sd %>% 
    		group_by(bin, !!sym(y_var), network, task) %>% 
    		summarize(cors = mean(cors)) %>% 
    		mutate(network = factor(network, levels = networks),
    		       cors = fisherz(cors))
	}
	
	if (!is.na(write_csv)) {
	    write.csv(sd, path(fig_save_root, write_csv, ext='csv'), row.names=FALSE)
	}
	
	# Determine scales
	if (!is.na(scales)) {
	    small <- floor(scales[1] * 100) / 100
	    big <- ceiling(scales[2] * 100) / 100
	} else {
    	small <- floor(min(pd$cors) * 100) / 100
    	big <- ceiling(max(pd$cors) * 100) / 100
	}
	
	# Whether to label the middle of scale
    if (label_middle) {
		sfg <- scale_fill_gradientn(colors = colors,
							 values = rescale(c(small, 0, big)),
							 breaks = c(small, 0, big),
							 limits = c(small, big))  
    } else {
		sfg <- scale_fill_gradientn(colors = colors,
							 values = rescale(c(small, 0, big)),
							 breaks = c(small, big),
							 limits = c(small, big))  
    }

	# Generate plot
	p <- pd %>%     
		ggplot(aes(x = bin, y = !!sym(y_var))) + 
		geom_tile(aes(fill = cors)) + 
		geom_point(data=ps, aes(x = bin, y = !!sym(y_var)), shape = 8, color = 'gold', size = 3) + 
		labs(
			x = x_label,
			y = y_label,
			fill = latex2exp::TeX('$z\\,[\\rho\\,]_{~~EEG, fMRI}$')
		) + 
	    sfg + 
		theme_bw() + 
		theme(strip.background = element_rect(fill = NA),
			  axis.ticks = element_blank(),
			  panel.grid = element_blank(),
			  legend.position = 'bottom',
			  text = element_text(size = overall_text),
			  axis.text = element_text(size = axis_text),
			  axis.text.x = element_text(angle = 45, hjust=1),
			  legend.text = element_text(size = legend_text, angle = 45, hjust=1),
			  legend.title = element_text(margin = margin(r = 30)))
		
	if (bytask) {
	    p <- p + facet_grid(task~network)
	} else {
	    p <- p + facet_wrap(~network)
	}
	if (!is.na(nrow)) p <- p + facet_wrap(~network, nrow = nrow)
	if (y_var == 'lag') p <- p + scale_y_continuous(breaks = seq(0, max(d$lag), 2), labels = seq(0, max(d$lag), 2)) 
	if (title == '') p <- p + theme(plot.title = element_blank())
	
	return (p)
}
