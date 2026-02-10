library(future.apply)
library(data.table)
library(tidyverse)
library(glue)
library(here)
library(fs)
library(arrow)

# This needs thorough commenting
# Skipping subject 23 session 2 for now
# (need to go back to the bids data and look at events.tsv)

# ---- GENERAL CORRELATION COMPUTATION FUNCTIONS ---- #

# --- FUNCTION DEFINITIONS --- #
# ---------------------------- #

# Calculate correlations

compute_correlations <- function(d, preprocess_root, by_task=FALSE) {
    
    setDT(d)
    
    # Make / refresh cache dir
    if (dir.exists(path(preprocess_root, 'cache'))) {
        dir_delete(path(preprocess_root, 'cache'))
    }
    
    dir.create(path(preprocess_root, 'cache'))
    
    # Identify EEG and fMRI columns
    first_eeg_col <- colnames(d)[which(grepl('^Fp1', colnames(d)))[1]]
    if (is.na(first_eeg_col)) stop('No EEG column starting with "Fp1" found.')
    eeg_cols <- colnames(d)[which(colnames(d) == first_eeg_col):(ncol(d))]
    fmri_cols <- colnames(d)[(which(colnames(d) == 'tr')+1):(which(colnames(d) == first_eeg_col)-1)]
    
    # Get networks from d
    # Make running by condition conditional on user input
    if (by_task) {
        condition_codes <- get_condition_codes(data_root, network=fmri_cols[1])
        d <- d[condition_codes, on = c('subject', 'session', 'run', 'tr'), nomatch = 0]
        # Rearrange columns
        d <- d[, c(colnames(condition_codes), colnames(d)[!colnames(d) %in% colnames(condition_codes)]), with = FALSE]
    }
    
    # Split data by subject, session, run
    if (by_task) {
        d[, run_id := paste(subject, session, run, condition, sep = '_')]
    } else {
        d[, run_id := paste(subject, session, run, sep = '_')]
    }
    
    # Get unique run_ids to distribute to workers
    run_ids <- unique(d$run_id)
    
    # Plan parallelization
    workers <- max(1, parallel::detectCores()-10)
    if (interactive() && .Platform$OS.type == 'unix') {
        # Interactive RStudio, use safer multisession
        plan(multisession, workers = 1)
    } else if (.Platform$OS.type == 'unix') {
        plan(multicore, workers=workers)
    } else {
        plan(multisession, workers = min(workers, 2))
    }
    
    # Increase globals size to avoid serialization errors
    options(future.globals.maxSize = 16 * 1024^3) # 16 GB
    
    print('Computing correlations...')
    
    future_lapply(run_ids, function(rid) {
        run_df <- d[run_id == rid]
        eeg_mat <- as.matrix(run_df[, ..eeg_cols])
        fmri_mat <- as.matrix(run_df[, ..fmri_cols])
        cor_mat <- cor(eeg_mat, fmri_mat, method = 'spearman')
        formatted <- format_cors(cor_mat, rid, by_task = by_task)
        write_feather(formatted, path(preprocess_root, glue('cache/{rid}.feather')))
        'Worker completed successfully'
    }, future.seed = FALSE)
    
}

format_cors <- function(data, label, by_task=FALSE) {
    
    # Parse metadata
    info <- str_split(label, '_')[[1]]
    subject <- info[1]
    session <- info[2]
    run <- info[3]
    if (by_task) task <- info[4]
    
    # Concatenate
    header <- data.frame(subject=subject, session=session, run=run, eeg_feature=rownames(data))
    if (by_task) header$task <- task
    out <- cbind(header, data)
    return(out)
}

correlation_pipeline <- function(d, script_root, ch_names, by_task=FALSE) {
    # d: merged_data.feather
    # script_root: The root of the preprocessing directory
    # ch_names: A string vector of channel names from MNE
    # by_task: Whether to compute correlations by task
    
    compute_correlations(d, script_root, by_task)
    
    # -- FORMAT RESULT -- #
    # For each (fmri_networks, eeg_features) matrix, transpose it, add it subject/session/run
    # info, and concatenate
    print('Formatting result...')
    rm(d)
    gc()
    cache_dir <- path(script_root, 'cache')
    files <- dir_ls(cache_dir, glob = '*.feather')
    
    d <- do.call(rbind, lapply(files, read_feather))
    
    d <- d %>% 
        separate(eeg_feature, into = c('channel', 'frequency', 'lag'), sep = '_') %>%
        mutate(frequency = as.numeric(frequency),
               lag = as.numeric(lag),
               channel = factor(channel, levels = ch_names))
    
    fmri_cols <- colnames(d)[(which(colnames(d) == 'lag')+1):(ncol(d))]
    fmri_cols <- fmri_cols[fmri_cols != 'task']
    
    d <- d %>% 
        gather(network, cors, all_of(fmri_cols)) 
    
    return(d)
    
}













# ---- FOR PARSING --BYTASK=true ---- #


parse_run <- function(combo, data_root) {
    subject <- combo$subject
    session <- combo$session
    network <- combo$network
    dir_path <- file.path(data_root, subject, session, 'func/conditions')
    
    # Import gradcpt
    pattern <- glue('{network}_.*{subject}_bld001.*')
    file <- list.files(path(dir_path, '../general'), pattern = pattern, full.names = TRUE)
    fmri <- sapply(readLines(file), as.numeric)
    gradcpt <- data.frame(run = 'run-001', tr = 1:length(fmri), condition = 'gradcpt', fmri = fmri)
    
    task_dfs <- list()
    # Import experience sampling / rest
    for (task in c('TP', 'Rest')) {
    
        pattern <- glue('{task}_{network}.*{subject}.*')
        files <- list.files(dir_path, pattern=pattern, full.names = TRUE)
        
        task_dfs[[task]] <- do.call(rbind, lapply(files, FUN = function(x) {
            # Make run string with three zero pads
            run_int <- as.integer(str_extract(x, 'bld(\\d+)', group = 1))
            run <- sapply(run_int, FUN = function(x) glue('run-00{x}'))
            fmri = sapply(readLines(x), as.numeric)
            d <- data.frame(run = run, tr = 1:length(fmri)) 
            d[[glue('{task}_fmri')]] <- fmri
            d
        }))
    }
    
    
    # -- Format task dfs -- 
    
    # Gradcpt
    gradcpt <- gradcpt %>% 
        mutate(condition =  ifelse(is.nan(fmri), 'head_motion', 'gradcpt')) %>% 
        select(-fmri) 
    
    # Experience sampling
    es <- cbind(task_dfs[[1]], rest_fmri = task_dfs[[2]][,'Rest_fmri'])
    
    es <- es %>% 
        mutate(condition = case_when(
            is.nan(TP_fmri) & is.nan(rest_fmri) ~ 'head_motion',
            is.nan(TP_fmri) & !is.nan(rest_fmri) ~ 'rest',
            !is.nan(TP_fmri) & is.nan(rest_fmri) ~ 'tp',
            .default = '-99'
        )) %>% 
        select(-TP_fmri, -rest_fmri) 
    
    out <- rbind(gradcpt, es)
    
    # Drop first lags
    num_lags <- suppressWarnings(as.integer(readLines(path(data_root, '../formatted/num_lags.txt'))))
    out <- out[num_lags:nrow(out),]
    
    # Drop head motion
    out <- out[out$condition != 'head_motion',]
    
    stem <- data.frame(subject=subject, session=session)    
    return(cbind(stem, out))
    
}

get_condition_codes <- function(data_root, network) {
    
    subjects <- list.dirs(path(data_root), recursive=FALSE, full.names = FALSE)
    sessions <- c('ses-001', 'ses-002')
    
    combos <- expand.grid(subject=subjects, session=sessions,  network=network)
    # Yank sub 23 ses 002 for now
    #combos <- combos[!(combos$subject == 'sub-023' & combos$session == 'ses-002'),]
    
    d <- do.call(rbind, lapply(1:nrow(combos), FUN=function(x) parse_run(combos[x, ], data_root)))
    
    return(d)
}


