% Analyze relationship between GradCPT RT variability and EEG signatures
% for DN and DAN

sublist = {'sub-002';'sub-003';'sub-004';'sub-006';'sub-007';'sub-008';'sub-009';'sub-010';...
    'sub-011';'sub-012';'sub-013';'sub-014';'sub-015';'sub-016';'sub-017';'sub-018';'sub-019';...
    'sub-020';'sub-021';'sub-023';'sub-024';'sub-025'};
sessions = {'ses-001';'ses-002'};
skip_sess1 = {'sub-003';'sub-025'};
skip_sess2 = {'sub-017'};

%% Set paths
bids_path = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/bids_data'];
eeg_signatures = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/correlation_data/eeg_signatures.csv'];
fmri_ts_path = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/MSHBM_sdc'];
eeg_signatures = importdata(eeg_signatures);
eeg_signature_headers = eeg_signatures.textdata(1,:);
eeg_signatures.textdata(1,:)=[];

%% Set params
TR = 2; % TR in seconds
nTR = 255;

%% Loop through subjects and sessions
% set up blank outputs
behav_onset = NaN(length(sublist),length(sessions));
OE_rate = NaN(length(sublist),length(sessions));
CE_rate = NaN(length(sublist),length(sessions));
n_BOLD_TRs = NaN(length(sublist),length(sessions));
n_EEG_samples = NaN(length(sublist),length(sessions));
n_RT_missing = NaN(length(sublist),length(sessions));
VTC_DNa_rho_corr = NaN(length(sublist),length(sessions));
VTC_DNa_dATNa_diff_corr = NaN(length(sublist),length(sessions));
VTC_DNa_dot_corr = NaN(length(sublist),length(sessions));
in_zone_DNa_rho = NaN(length(sublist),length(sessions));
out_zone_DNa_rho = NaN(length(sublist),length(sessions));
in_zone_DNa_dATNa_diff = NaN(length(sublist),length(sessions));
out_zone_DNa_dATNa_diff = NaN(length(sublist),length(sessions));
DNa_BOLD_EEG_rho_corr = NaN(length(sublist),length(sessions));
DNa_dATNa_rho_corr = NaN(length(sublist),length(sessions));
DNa_dATNa_dot_corr = NaN(length(sublist),length(sessions));
DNa_BOLD_EEG_dot_corr = NaN(length(sublist),length(sessions));
DNa_dATNa_BOLD_corr = NaN(length(sublist),length(sessions));
VTC_DNa_BOLD_corr = NaN(length(sublist),length(sessions));
mean_DNa_BOLD_before_mt = NaN(length(sublist),length(sessions));
mean_DNa_BOLD_after_mt = NaN(length(sublist),length(sessions));
mean_DNa_EEG_before_mt = NaN(length(sublist),length(sessions));
mean_DNa_EEG_after_mt = NaN(length(sublist),length(sessions));
DNa_BOLD_after_minus_before_mt = NaN(length(sublist),length(sessions));
DNa_EEG_after_minus_before_mt = NaN(length(sublist),length(sessions));
VTC_dATNa_rho_corr = NaN(length(sublist),length(sessions));
VTC_dATNa_dot_corr = NaN(length(sublist),length(sessions));
dATNa_BOLD_EEG_rho_corr = NaN(length(sublist),length(sessions));
dATNa_BOLD_EEG_dot_corr = NaN(length(sublist),length(sessions));
VTC_dATNa_BOLD_corr = NaN(length(sublist),length(sessions));
dATNa_BOLD_after_minus_before_mt = NaN(length(sublist),length(sessions));
dATNa_EEG_after_minus_before_mt = NaN(length(sublist),length(sessions));

VTC_FPCNb_rho_corr = NaN(length(sublist),length(sessions));
VTC_FPCNb_dot_corr = NaN(length(sublist),length(sessions));
FPCNb_BOLD_EEG_rho_corr = NaN(length(sublist),length(sessions));
FPCNb_BOLD_EEG_dot_corr = NaN(length(sublist),length(sessions));
VTC_FPCNb_BOLD_corr = NaN(length(sublist),length(sessions));
FPCNb_BOLD_after_minus_before_mt = NaN(length(sublist),length(sessions));
FPCNb_EEG_after_minus_before_mt = NaN(length(sublist),length(sessions));

% loop through
for i = 1:length(sublist)
    for j = 1:length(sessions)
        if ~(sum(contains(skip_sess2,sublist{i}))>0 && j==2) && ~(sum(contains(skip_sess1,sublist{i}))>0 && j==1)
        %% Compute the VTC and assign to each TR
        % Load GradCPT behav data
        behav = [];
        behav = dir([bids_path filesep sublist{i} filesep sessions{j} '/func/sub*GradCPT*tsv']);
        behav = importdata([behav.folder filesep behav.name]);
        behav(1,:) = []; % delete first trial (no response)
        behav(end,:) = []; % delete last trial (no response)
        behav_onset(i,j) = behav(1,1);
        % Load BOLD time series
        DNa_BOLD = []; dATNa_BOLD = []; FPCNb = [];
        DNa_BOLD = load([fmri_ts_path '/' sublist{i} '/' sessions{j} '/MSHBM_ts/DNa_ts_MSHBM_' sublist{i} '_bld001_highpass.txt']);
        dATNa_BOLD = load([fmri_ts_path '/' sublist{i} '/' sessions{j} '/MSHBM_ts/dATNa_ts_MSHBM_' sublist{i} '_bld001_highpass.txt']);
        FPCNb_BOLD = load([fmri_ts_path '/' sublist{i} '/' sessions{j} '/MSHBM_ts/FPCNb_ts_MSHBM_' sublist{i} '_bld001_highpass.txt']);
        % Get OE and CE rates
        OE_inds = []; CE_inds = []; CO_inds =[];
        OE_inds = find(behav(:,3)==2 & behav(:,4)==0);
        OE_rate(i,j) = length(OE_inds)/length(find(behav(:,3)==2));
        CE_inds = find(behav(:,3)==1 & behav(:,4)~=0);
        CE_rate(i,j) = length(CE_inds)/length(find(behav(:,3)==1));
        % Covert OEs and COs to NaN for RT
        RTs_all = [];
        CO_inds = find(behav(:,3)==1 & behav(:,4)==0);
        RTs_all = behav(:,7);
        RTs_all([CO_inds; OE_inds]) = NaN;
        % Compute VTC
        VTC = []; abs_VTC = [];
        meanRT = nanmean(RTs_all);
        stdRT = nanstd(RTs_all);
        VTC = (RTs_all-meanRT)/stdRT;
        abs_VTC = abs(VTC);
        % IGNORE MULTI-TRIAL WINDOWS WITH NO RESPONSES?
        % interpolate NaNs
        abs_VTC_interp = [];
        n_RT_missing(i,j) = sum(isnan(abs_VTC));
        abs_VTC_interp = inpaint_nans(abs_VTC);
        abs_VTC_interp = abs(abs_VTC_interp);
        % Get TR time course
        TRs = 0; TR_step = TR;
        for k=1:nTR
            TRs = [TRs; TR_step];
            TR_step = TR_step +TR;
        end  
        % get mean VTC within 2-second window prior to each TR onset
        RT_onsets = []; VTC_for_TR = [];
        RT_onsets = behav(:,1);
        for k=1:length(TRs)
            RT_inds_for_TR = find(TRs(k)-RT_onsets>0 & TRs(k)-RT_onsets<TR);
            if ~isempty(RT_inds_for_TR)
                VTC_for_TR(k) = mean(abs_VTC_interp(RT_inds_for_TR));
            else
                VTC_for_TR(k) = NaN;
            end
        end
        % get timing of mountain events
        mt_ind = find(behav(:,3)==1);
        mt_onsets_sec = behav(mt_ind,1);
        mt_onsets_TR = mt_onsets_sec/TR; % convert to TR time
        %% Get DNa and dATNa signature values per TR
        DNa_eeg_sig_inds = find(strcmp(eeg_signatures.textdata(:,1),sublist{i}) & strcmp(eeg_signatures.textdata(:,2),sessions{j}) & strcmp(eeg_signatures.textdata(:,4),'DNa'));
        dATNa_eeg_sig_inds = find(strcmp(eeg_signatures.textdata(:,1),sublist{i}) & strcmp(eeg_signatures.textdata(:,2),sessions{j}) & strcmp(eeg_signatures.textdata(:,4),'dATNa'));
        FPCNb_eeg_sig_inds = find(strcmp(eeg_signatures.textdata(:,1),sublist{i}) & strcmp(eeg_signatures.textdata(:,2),sessions{j}) & strcmp(eeg_signatures.textdata(:,4),'FPCNb'));
        TRs_included = str2double(eeg_signatures.textdata(DNa_eeg_sig_inds,3));
        DNa_TRs_sec = (str2double(eeg_signatures.textdata(DNa_eeg_sig_inds,3)))*TR;
        %check_DNa_TRs_sec{i,j} = DNa_TRs_sec; % CHECK
        DNa_EEG_rho = eeg_signatures.data(DNa_eeg_sig_inds,2);
        dATNa_EEG_rho = eeg_signatures.data(dATNa_eeg_sig_inds,2);
        FPCNb_EEG_rho = eeg_signatures.data(FPCNb_eeg_sig_inds,2);
        DNa_dATN_rho_diff = DNa_EEG_rho-dATNa_EEG_rho;
        DNa_EEG_dot = eeg_signatures.data(DNa_eeg_sig_inds,1);
        dATNa_EEG_dot = eeg_signatures.data(dATNa_eeg_sig_inds,1);
        FPCNb_EEG_dot = eeg_signatures.data(FPCNb_eeg_sig_inds,1);
        VTC_val = [];
        for k=1:length(DNa_TRs_sec)
            VTC_ind = find(DNa_TRs_sec(k)==TRs);
            VTC_val(k) = VTC_for_TR(VTC_ind);
        end
        % Correlate VTC vs EEG DNa and dATNa
        VTC_DNa_rho_corr(i,j) = corr(DNa_EEG_rho,VTC_val','rows','complete');
        VTC_dATNa_rho_corr(i,j) = corr(dATNa_EEG_rho,VTC_val','rows','complete');
        VTC_FPCNb_rho_corr(i,j) = corr(dATNa_EEG_rho,VTC_val','rows','complete');
        VTC_DNa_dot_corr(i,j) = corr(DNa_EEG_dot,VTC_val','rows','complete');
        VTC_dATNa_dot_corr(i,j) = corr(dATNa_EEG_dot,VTC_val','rows','complete');
        VTC_FPCNb_dot_corr(i,j) = corr(dATNa_EEG_dot,VTC_val','rows','complete');
        VTC_DNa_dATNa_diff_corr(i,j) = corr(DNa_dATN_rho_diff,VTC_val','rows','complete');
        % Correlate VTC vs BOLD DNa and dATNa
        DNa_BOLD(1:10) = []; % delete first 10 TRs (to match EEG data)
        dATNa_BOLD(1:10) = [];
        FPCNb_BOLD(1:10) = [];
        DNa_BOLD(isnan(DNa_BOLD))=[]; % delete BOLD TRs with head motion
        dATNa_BOLD(isnan(dATNa_BOLD))=[];
        FPCNb_BOLD(isnan(FPCNb_BOLD))=[];
        DNa_BOLD_HRF = DNa_BOLD;
        DNa_BOLD_HRF(1:3) = []; % Shift BOLD time series by 6 seconds to account for HRF
        DNa_BOLD_HRF = [DNa_BOLD_HRF; NaN; NaN; NaN];
        dATNa_BOLD_HRF = dATNa_BOLD;
        dATNa_BOLD_HRF(1:3) = []; % Shift BOLD time series by 6 seconds to account for HRF
        dATNa_BOLD_HRF = [dATNa_BOLD_HRF; NaN; NaN; NaN];
        FPCNb_BOLD_HRF = FPCNb_BOLD;
        FPCNb_BOLD_HRF(1:3) = []; % Shift BOLD time series by 6 seconds to account for HRF
        FPCNb_BOLD_HRF = [FPCNb_BOLD_HRF; NaN; NaN; NaN];
        n_BOLD_TRs(i,j) = length(DNa_BOLD_HRF);
        n_EEG_samples(i,j) = length(DNa_EEG_rho);
        if length(DNa_BOLD_HRF)>length(DNa_EEG_rho) % account for runs with missing EEG TRs
           TR_diff = n_BOLD_TRs(i,j)-n_EEG_samples(i,j);
           DNa_BOLD_HRF(end-TR_diff+1:end) = []; 
           dATNa_BOLD_HRF(end-TR_diff+1:end) = []; 
           FPCNb_BOLD_HRF(end-TR_diff+1:end) = []; 
           n_BOLD_TRs(i,j) = length(DNa_BOLD_HRF);
        end
        VTC_DNa_BOLD_corr(i,j) = corr(DNa_BOLD_HRF,VTC_val','rows','complete');
        VTC_dATNa_BOLD_corr(i,j) = corr(dATNa_BOLD_HRF,VTC_val','rows','complete');
        VTC_FPCNb_BOLD_corr(i,j) = corr(FPCNb_BOLD_HRF,VTC_val','rows','complete');
        % Correlate DNa EEG vs dATNa EEG
        DNa_dATNa_rho_corr(i,j) = corr(DNa_EEG_rho,dATNa_EEG_rho);
        DNa_dATNa_dot_corr(i,j) = corr(DNa_EEG_dot,dATNa_EEG_dot);
        % Correlate DNa BOLD vs dATNa BOLD
        DNa_dATNa_BOLD_corr(i,j) = corr(DNa_BOLD,dATNa_BOLD,'rows','complete');
        % Median split (or percentile cutoff) on VTC
        high_cutoff = prctile(VTC_val,75);
        low_cutoff = prctile(VTC_val,25);
        in_zone_inds = find(VTC_val<low_cutoff);
        out_zone_inds = find(VTC_val>high_cutoff);
        in_zone_DNa_rho(i,j) = mean(DNa_EEG_rho(in_zone_inds));
        out_zone_DNa_rho(i,j) = mean(DNa_EEG_rho(out_zone_inds));
        in_zone_dATNa_rho(i,j) = mean(dATNa_EEG_rho(in_zone_inds));
        out_zone_dATNa_rho(i,j) = mean(dATNa_EEG_rho(out_zone_inds));
        in_zone_DNa_dATNa_diff(i,j) = mean(DNa_dATN_rho_diff(in_zone_inds));
        out_zone_DNa_dATN_diff(i,j) = mean(DNa_dATN_rho_diff(out_zone_inds));
        % Correlate BOLD vs EEG signature
        DNa_BOLD_EEG_rho_corr(i,j) = corr(DNa_BOLD_HRF,DNa_EEG_rho,'rows','complete');
        dATNa_BOLD_EEG_rho_corr(i,j) = corr(dATNa_BOLD_HRF,dATNa_EEG_rho,'rows','complete');
        FPCNb_BOLD_EEG_rho_corr(i,j) = corr(FPCNb_BOLD_HRF,dATNa_EEG_rho,'rows','complete');
        DNa_BOLD_EEG_dot_corr(i,j) = corr(DNa_BOLD_HRF,DNa_EEG_dot,'rows','complete');
        dATNa_BOLD_EEG_dot_corr(i,j) = corr(dATNa_BOLD_HRF,dATNa_EEG_dot,'rows','complete');
        FPCNb_BOLD_EEG_dot_corr(i,j) = corr(FPCNb_BOLD_HRF,dATNa_EEG_dot,'rows','complete');
        % Compare BOLD before vs. after mountain onsets
        mt_onsets_TR = unique(round(mt_onsets_TR)); % identify TRs adjacent to mt onsets
        mt_rm_ind = find(diff(mt_onsets_TR)==1); % remove mountain TRs preceded by prior mountain TRs
        mt_onsets_TR(mt_rm_ind+1)=[];
        DNa_BOLD_before_mt = []; DNa_BOLD_after_mt = [];
        DNa_EEG_before_mt = []; DNa_EEG_after_mt = [];
        dATNa_BOLD_before_mt = []; dATNa_BOLD_after_mt = [];
        dATNa_EEG_before_mt = []; dATNa_EEG_after_mt = [];
        FPCNb_BOLD_before_mt = []; FPCNb_BOLD_after_mt = [];
        FPCNb_EEG_before_mt = []; FPCNb_EEG_after_mt = [];
        DNa_BOLD_z = zscore(DNa_BOLD); % z score BOLD time series
        dATNa_BOLD_z = zscore(dATNa_BOLD);
        FPCNb_BOLD_z = zscore(FPCNb_BOLD);
        for k = 1:length(mt_onsets_TR)
            before_mt_TR_ind =  [find(TRs_included==mt_onsets_TR(k))-1]; % 1 TR before
            after_mt_TR_inds_EEG = [find(TRs_included==mt_onsets_TR(k))+0 find(TRs_included==mt_onsets_TR(k))+1]; % 0-1 TRs after
            after_mt_TR_inds_BOLD = [find(TRs_included==mt_onsets_TR(k))+3 find(TRs_included==mt_onsets_TR(k))+4]; % 3-4  TRs after for HRF
            if ~isempty(before_mt_TR_ind) && mt_onsets_TR(k)>11 && isempty(find(before_mt_TR_ind==0))==1 && sum(after_mt_TR_inds_BOLD>length(DNa_BOLD))==0
               % ignore mountains occuring at first TR of task or
               % surrounding missing TRs due to censoring
               DNa_BOLD_before_mt = [DNa_BOLD_before_mt; mean(DNa_BOLD_z(before_mt_TR_ind))];
               DNa_BOLD_after_mt = [DNa_BOLD_after_mt; mean(DNa_BOLD_z(after_mt_TR_inds_BOLD))];
               DNa_EEG_before_mt = [DNa_EEG_before_mt; mean(DNa_EEG_rho(before_mt_TR_ind))];
               DNa_EEG_after_mt = [DNa_EEG_after_mt; mean(DNa_EEG_rho(after_mt_TR_inds_EEG))];
               dATNa_BOLD_before_mt = [dATNa_BOLD_before_mt; mean(dATNa_BOLD_z(before_mt_TR_ind))];
               dATNa_BOLD_after_mt = [dATNa_BOLD_after_mt; mean(dATNa_BOLD_z(after_mt_TR_inds_BOLD))];
               dATNa_EEG_before_mt = [dATNa_EEG_before_mt; mean(dATNa_EEG_rho(before_mt_TR_ind))];
               dATNa_EEG_after_mt = [dATNa_EEG_after_mt; mean(dATNa_EEG_rho(after_mt_TR_inds_EEG))];
               FPCNb_BOLD_before_mt = [FPCNb_BOLD_before_mt; mean(FPCNb_BOLD_z(before_mt_TR_ind))];
               FPCNb_BOLD_after_mt = [FPCNb_BOLD_after_mt; mean(FPCNb_BOLD_z(after_mt_TR_inds_BOLD))];
               FPCNb_EEG_before_mt = [FPCNb_EEG_before_mt; mean(FPCNb_EEG_rho(before_mt_TR_ind))];
               FPCNb_EEG_after_mt = [FPCNb_EEG_after_mt; mean(FPCNb_EEG_rho(after_mt_TR_inds_EEG))];
            end
        end
        % Compare EEG signature before vs. after mountain onsets
        DNa_BOLD_after_minus_before_mt(i,j) = mean(DNa_BOLD_after_mt-DNa_BOLD_before_mt);
        DNa_EEG_after_minus_before_mt(i,j) = mean(DNa_EEG_after_mt-DNa_EEG_before_mt);
        dATNa_BOLD_after_minus_before_mt(i,j) = mean(dATNa_BOLD_after_mt-dATNa_BOLD_before_mt);
        dATNa_EEG_after_minus_before_mt(i,j) = mean(dATNa_EEG_after_mt-dATNa_EEG_before_mt);
        FPCNb_BOLD_after_minus_before_mt(i,j) = mean(FPCNb_BOLD_after_mt-FPCNb_BOLD_before_mt);
        FPCNb_EEG_after_minus_before_mt(i,j) = mean(FPCNb_EEG_after_mt-FPCNb_EEG_before_mt);
        mean_DNa_BOLD_before_mt(i,j) = mean(DNa_BOLD_before_mt);
        mean_DNa_BOLD_after_mt(i,j) = mean(DNa_BOLD_after_mt);        
        end
    end
end

% compute mean corrs per subject across session
mean_VTC_DNa_rho_corr = nanmean(VTC_DNa_rho_corr,2);
mean_VTC_dATNa_rho_corr = nanmean(VTC_dATNa_rho_corr,2);
mean_VTC_DNa_dATN_diff_corr = nanmean(VTC_DNa_dATNa_diff_corr,2);

%table(sublist,behav_onset)