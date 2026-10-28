%% Get correlations of each freq and lag (channel-averaged) with BOLD 

% setup
networks = {'DNa'; 'DNb'};
data_dir = '/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/formatted/full';
roi_dir = '/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/preprocessed';
sublist = dir([data_dir '/sub*.mat']);
sublist = {sublist.name};
session_names = {'ses-001';'ses-002'};
run_names = {'bld001'; 'bld002'; 'bld003'; 'bld004'};
lags = 11;
freqs = 40;
n_lag_freq_pairs = lags*freqs;

for n=1:length(networks)
    ROI = networks{n};
% loop through subjects
r_sub_avg = [];
for i=1:length(sublist)
    data=[];
    load([data_dir filesep sublist{i}]);
    sessions = fieldnames(data);
    % loop through sessions
    r_sess_avg = [];
    r_allsessions = [];
    for j=1:length(sessions)
        session_data = [];
        session_data = data.(sessions{j});
        runs = fieldnames(session_data);
        % loop through runs
        r_avg_allruns = [];
        for k=1:length(runs)
            run_data = []; BOLD_ts = []; corr_eeg_fmri = [];
            run_data = session_data.(runs{k});
            % compute correlations between BOLD and all channels, frequencies and lags
            if i==3 && j==1 % account for skipped first run in sub-003 ses-001
                 BOLD_ts = load([roi_dir '/' sublist{i}(1:7) '/' session_names{j} '/func/general/'...
                    ROI '_ts_MSHBM_' sublist{i}(1:7) '_' run_names{k+1} '_highpass.txt']);
            else
               BOLD_ts = load([roi_dir '/' sublist{i}(1:7) '/' session_names{j} '/func/general/'...
                    ROI '_ts_MSHBM_' sublist{i}(1:7) '_' run_names{k} '_highpass.txt']);
            end
            
            BOLD_ts(1:10) = []; % remove first 10 TRs
            BOLD_ts(isnan(BOLD_ts))=[]; % remove NaNs (censored)
            if length(BOLD_ts)>size(run_data.X,1) % trim BOLD if needed
                BOLD_eeg_diff=length(BOLD_ts)-size(run_data.X,1);
                BOLD_ts(end-BOLD_eeg_diff+1:end)=[];
            end
            %BOLD_ts = run_data.y.(ROI); % CHANGE THIS FOR NULL MODELS
            
            corr_eeg_fmri = corr(BOLD_ts,run_data.X,'type','Spearman');
            % average across channels for each frequency-lag pair
            r_chan_avg = [];
            for l=1:n_lag_freq_pairs
                curr_idx = [];
                curr_idx = l:n_lag_freq_pairs:size(run_data.X,2);
                r_chan_avg(l) = mean(corr_eeg_fmri(curr_idx));                
            end
            % average within freqency bands
            r_freq_avg = []; % for storing freq averaged values
            bands = {[1:3]; [4:7]; [8:12]; [13:20]; [21:30]; [31:40]};
            for l=1:length(bands)
                curr_band = []; curr_vals = [];
                curr_band = bands{l}; % extract all values for this band
                if l>1 
                    curr_vals = r_chan_avg((curr_band(1)*lags+1):(curr_band(end)*lags));
                else
                    curr_vals = r_chan_avg(1:curr_band(end)*lags);
                end
            curr_matrix = reshape(curr_vals, lags, []); % average for every lag for this band
            r_freq_avg = [r_freq_avg; mean(curr_matrix,2)];
            end
            % average across lags 0-10 seconds for each frequency
            r_lag_avg = [];
            freq_start_idx = 1;
            for l=1:length(bands)
                r_lag_avg = [r_lag_avg; mean(r_freq_avg(freq_start_idx:freq_start_idx+5))];
                freq_start_idx = freq_start_idx + lags; 
            end 
            r_avg_allruns = [r_avg_allruns r_lag_avg];
            display(['done ' sublist{i} ' ' sessions{j} ' ' runs{k}]);
        end
        % average across runs within session
        r_sess_avg = mean(r_avg_allruns,2);
        r_allsessions = [r_allsessions r_sess_avg];
        %r_sess_avg = r_lag_avg;
        %r_sess_avg = squeeze(mean(r_freq_avg,2));
    end
    % average across sessions within subject
    r_sub_avg = [r_sub_avg mean(r_allsessions,2)];
    r_sub_avg_allnetworks{n} = r_sub_avg;
    %r_sub_avg = [r_sub_avg; mean(r_sess_avg)];
end
end

% get differences between networks within subjects (for null model)
network_diff = r_sub_avg_allnetworks{1,1}-r_sub_avg_allnetworks{1,2};
abs_network_diff = abs(r_sub_avg_allnetworks{1,1}-r_sub_avg_allnetworks{1,2});
mkdir('results');
save(['results/real_distribution_signed.mat'],'network_diff');
save(['results/real_distribution.mat'],'abs_network_diff');


