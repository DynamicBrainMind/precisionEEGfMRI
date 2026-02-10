
plot_by_task <- function(d_task, networks, bands, axis_text=14, overall_text=18,
                         scales=NA, title='', x_label=NA, y_label=NA, by_channels=FALSE,
                         colors=NA, nlag_s=NA, write_csv=NA) {
    
    # Configure inputs
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
    
    if (!is.na(nlag_s)) d_task <- d_task[d_task$lag <= nlag_s,]
    
    # Summarize data
    sd <- d_task %>% 
        group_by(subject, network, frequency, !!sym(y_var), task) %>% 
        summarize(cors = mean(cors)) 
    pd <- sd %>% 
        group_by(network, frequency, !!sym(y_var), task) %>% 
        summarize(cors = mean(cors)) 
    
    write.csv(sd, path(fig_save_root, write_csv, ext='.csv'))
    
    # Set scale constraints
    if (!is.na(scales)) {
        small <- floor(scales[1] * 100) / 100
        big <- ceiling(scales[2] * 100) / 100
    } else {
        small <- floor(min(pd$cors)*100)/100
        big <- ceiling(max(pd$cors)*100)/100
    }
    
    p <- pd %>% 
        filter(network %in% !!networks) %>% 
        ggplot(aes(x = frequency, y = !!sym(y_var))) + 
        geom_tile(aes(fill = cors)) + 
        facet_grid(task ~ network) + 
        labs(
            title = title,
            x = x_label,
            y = y_label,
            fill = latex2exp::TeX('$\\rho_{~~EEG,fMRI}$')
        ) + 
        scale_fill_gradientn(colors = colors,
                             values = rescale(c(small, 0, big)),
                             limits = c(small, big),
                             breaks = c(small, 0, big),
                             labels = c(small, 0, big)) + 
        theme_bw() + 
        theme(panel.grid = element_blank(),
              axis.ticks = element_blank(),
              strip.background = element_rect(fill = NA, color = 'black'),
              legend.title = element_text(margin = margin(r = 30)),
              text = element_text(size = overall_text),
              axis.text = element_text(size = axis_text),
              legend.position = 'bottom')
    
    if(maybe_log_scale(d$frequency)) p <- p + scale_x_log10()
    if (y_var == 'lag') p <- p + scale_y_continuous(breaks = seq(0, max(d$lag), 2), labels = seq(0, max(d$lag), 2)) 
    if (title == '') p <- p + theme(plot.title = element_blank())
    
    return(p)
    
}
