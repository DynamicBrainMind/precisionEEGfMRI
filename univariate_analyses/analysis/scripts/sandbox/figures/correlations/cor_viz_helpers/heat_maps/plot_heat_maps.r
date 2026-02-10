maybe_log_scale <- function(x, tol = 1e-6) {
    # Check whether a vector x is log scaled
    x <- unique(x)
    n <- length(x)
    ideal <- 10^seq(log10(min(x)), log10(max(x)), length.out = n)
    return(all(abs(sort(x) - ideal) < tol))
}

plot_heat <- function(d, networks, scales = NA,
                      overall_text=18, axis_text=16, x_label=NA, y_label=NA,
                      title='', by_channels=FALSE, colors=NA, nlags_s=NA, nrow=NA,
                      write_csv=NA) {
    #' Plot Heatmap of Mean Correlations by Frequency and Lag
    #'
    #' Generates a faceted heatmap of mean correlation values (`cors`) across
    #' subjects for selected networks, as a function of frequency and lag.
    #'
    #' @param d A data frame containing columns `subject`, `lag`, `frequency`, `network`, and `cors`.
    #' @param networks A character vector of network names to include in the plot.
    #' @param scales Optional numeric vector `c(min, max)` to set heatmap fill limits. Defaults to automatic scaling.
    #' @param overall_text Numeric, base font size for overall plot text. Default is 18.
    #' @param axis_text Numeric, font size for axis labels. Default is 16.
    #'
    #' @return A `ggplot` object representing the faceted heatmap.
    #'
    #' @details
    #' - Frequencies are binned and labeled according to standard EEG bands:
    #'   Delta, Theta, Alpha, Beta, Gamma.
    #' - Mean correlations are computed first per subject, then averaged across subjects.
    #' - Heatmap fill uses a diverging red-blue color scale centered at 0.
    #' - Optionally applies a log scale to the x-axis if `maybe_log_scale()` returns TRUE.
    #'
    #' @examples
    #' plot_heat(df, networks = c("DMN", "FPN"), scales = c(-0.2, 0.3))
    
    y_var <- ifelse(by_channels, 'channel', 'lag')
    y_label <- ifelse(!is.na(x_label), x_label,
                      ifelse(y_var == 'channel', 'Channel', 'Lag relative to BOLD (s)'))
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
	bins <- levels(cut(d$frequency, breaks=breaks))
	labels <- paste(labels, bins, sep=' ')

	# Prep data
	if (!is.na(nlags_s)) d <- d[d$lag <= nlags_s,]
	
	sd <- d %>% 
		filter(network %in% networks) %>%
		group_by(subject, !!sym(y_var), frequency, network) %>% 
		summarize(cors = mean(cors)) 
    pd <- sd %>% 	
		group_by(!!sym(y_var), frequency, network) %>% 
		summarize(cors = mean(cors)) 
	
	if (!is.na(write_csv)) {
	    write.csv(sd, path(fig_save_root, write_csv, ext='csv'), row.names = FALSE) # Pulls fig_save_root from global...
	}

	# Set scale constraints
	if (!is.na(scales)) {
	    small <- floor(scales[1] * 100) / 100
	    big <- ceiling(scales[2] * 100) / 100
	} else {
    	small <- floor(min(pd$cors)*100)/100
    	big <- ceiling(max(pd$cors)*100)/100
	}

	# Plot
	p1 <- pd %>% 
		mutate(network = factor(network, levels = networks)) %>%
		ggplot(aes(x = frequency, y = !!sym(y_var))) +
	    geom_tile(aes(fill = cors)) +
		facet_wrap(~network) +
		scale_fill_gradientn(colors = colors,
							 values = rescale(c(small, 0, big)),
							 limits = c(small, big),
							 breaks = c(small, 0, big),
							 labels = c(small, 0, big)) + 
		labs(
		    title = title,
			x = 'Frequency (Hz)',
			y = y_label,
			fill = latex2exp::TeX('$\\rho_{~~EEG, fMRI}$')
		) + 
		theme_bw() + 
		theme(strip.background = element_rect(fill = NA),
			  panel.grid = element_blank(),
			  axis.ticks = element_blank(),
			  axis.text = element_text(size = axis_text),
			  #legend.position = 'right',
			  legend.position = 'bottom',
			  legend.direction = 'horizontal',
			  text = element_text(size=overall_text))
	
	if (!is.na(nrow)) p1 <- p1 + facet_wrap(~network, nrow=nrow)
	if (y_var == 'lag') p1 <- p1 + scale_y_continuous(breaks = seq(0, max(d$lag), 2), labels = seq(0, max(d$lag), 2)) 
	if (maybe_log_scale(pd$frequency)) p1 <- p1 + scale_x_log10()
	if (title == '') p1 <- p1 + theme(plot.title = element_blank())

	return(p1)

}
