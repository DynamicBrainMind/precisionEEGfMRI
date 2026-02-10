% Use experience-sampling TR condition labels (from
% mk_ES_condition_labels.m) to make condition-specific versions of ROI time
% series
% also computes n TRs retained for analysis in each condition

%% settings
mshbm_dir = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/MSHBM_sdc'];
arealmshbm_dir = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/ArealMSHBM_sdc'];
sublist = dir('/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/CBIG_preproc_fmap_highpass/sub*');
sublist = {sublist.name};
sessions = {'ses-001'; 'ses-002'};
runs = {'bld002';'bld003';'bld004'};
nparcels='300';

for i = 1:length(sublist)
    sub = sublist{i};
    for j = 1:length(sessions)
        sess = sessions{j};
        for k = 1:length(runs)
        % Get condition TR label files
        Condition_file = dir([mshbm_dir '/' sub '/' sess '/MSHBM_ts/ConditionLabel*' runs{k} '.txt']);
        if ~isempty(Condition_file)
        cond_labels = load([Condition_file.folder filesep Condition_file.name]);
        non_rest_inds = find(cond_labels~=1);
        non_tp_inds = find(cond_labels~=2);
        % Get all MSHBM and ArealMSHBM ROI time series files for each ES run, label
        % conditions, and save outputs
        BOLD_ts_files = []; Areal_files = [];
        BOLD_ts_files = dir([mshbm_dir '/' sub '/' sess '/MSHBM_ts/*' runs{k} '_highpass*']);
        BOLD_ts_files(startsWith({BOLD_ts_files.name},'Rest_')) = []; % ignore previously created Rest files
        BOLD_ts_files(startsWith({BOLD_ts_files.name},'TP_')) = [];
        Areal_files = dir([arealmshbm_dir '/' sub '/' sess '/' nparcels '/ArealMSHBM_ts/*' runs{k} '*']);
        Areal_files(startsWith({Areal_files.name},'Rest_')) = [];
        Areal_files(startsWith({Areal_files.name},'TP_')) = [];
        for m = 1:length(BOLD_ts_files)
            BOLD_ts = []; rest_BOLD_ts = []; tp_BOLD_ts = [];
            BOLD_ts = load([BOLD_ts_files(1).folder filesep BOLD_ts_files(m).name]);
            rest_BOLD_ts = BOLD_ts; tp_BOLD_ts = BOLD_ts;
            rest_BOLD_ts(non_rest_inds) = NaN;
            tp_BOLD_ts(non_tp_inds) = NaN;
            dlmwrite([mshbm_dir '/' sub '/' sess '/MSHBM_ts/Rest_' BOLD_ts_files(m).name],rest_BOLD_ts);
            dlmwrite([mshbm_dir '/' sub '/' sess '/MSHBM_ts/TP_' BOLD_ts_files(m).name],tp_BOLD_ts);
        end
        for m = 1:length(Areal_files)-2 % remove last 2 files (all 300 regions)
            BOLD_ts = []; rest_BOLD_ts = []; tp_BOLD_ts = [];
            BOLD_ts = load([Areal_files(1).folder filesep Areal_files(m).name]);
            rest_BOLD_ts = BOLD_ts; tp_BOLD_ts = BOLD_ts;
            rest_BOLD_ts(non_rest_inds) = NaN;
            tp_BOLD_ts(non_tp_inds) = NaN;
            dlmwrite([arealmshbm_dir '/' sub '/' sess '/' nparcels '/ArealMSHBM_ts/Rest_' Areal_files(m).name],rest_BOLD_ts);
            dlmwrite([arealmshbm_dir '/' sub '/' sess '/' nparcels '/ArealMSHBM_ts/TP_' Areal_files(m).name],tp_BOLD_ts);
        end
        % count n total TRs retained in experience sampling runs (excluding head motion and ignored
        % samples)
        TRs_total{i}(j,k) = sum(~isnan(BOLD_ts));
        TRs_rest{i}(j,k) = sum(~isnan(rest_BOLD_ts));
        TRs_TP{i}(j,k) = sum(~isnan(tp_BOLD_ts));
        else
            display(['skipping ' sub ' ' sess ' ' runs{k} ' due to missing data'])
        end
        end
    end   
end

% compute total TRs per session
TRs_total_sess_ES = []; TRs_rest_sess_ES = []; TRs_TP_sess_ES = []; TRs_GradCPT = [];
TRs_all = NaN(length(sessions),length(sublist));
for i = 1:length(sublist)
    for j = 1:length(sessions)
        if size(TRs_total{i})>=j
        TRs_total_sess_ES = [TRs_total_sess_ES; sum(TRs_total{i}(j,:))];
        TRs_rest_sess_ES(i,j) = sum(TRs_rest{i}(j,:));
        %TRs_rest_sess_ES = [TRs_rest_sess_ES; sum(TRs_rest{i}(j,:))];
        TRs_TP_sess_ES(i,j) = sum(TRs_TP{i}(j,:));
        %TRs_TP_sess_ES = [TRs_TP_sess_ES; sum(TRs_TP{i}(j,:))];
        % load a sample GradCPT BOLD time series (with head motion
        % censored)
        BOLD_ts_files = dir([mshbm_dir '/' sublist{i} '/' sessions{j} '/MSHBM_ts/*bld001' '_highpass*']);
        BOLD_ts = load([BOLD_ts_files(1).folder filesep BOLD_ts_files(1).name]);
        TRs_GradCPT = sum(~isnan(BOLD_ts));
        TRs_all(j,i) = TRs_GradCPT + sum(TRs_total{i}(j,:));
        end
    end
end

