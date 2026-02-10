% For experience sampling data, label each TRs as rest (1) vs thought probe (2) and save to a vector

%% settings
bids_dir = '/Users/ak4379/Documents/data/R21_EEG-fMRI/bids_data';
mshbm_dir = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/MSHBM_sdc'];
%sublist = dir('/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/CBIG_preproc_fmap_highpass/sub*');
%sublist = {sublist.name};
sublist = {'sub-023'};
sessions = {'ses-001'; 'ses-002'};
mshbm_es_run_labels = {'bld002';'bld003';'bld004'};

%% Loop through subjects, sessions and runs: label TRs based on condition
 for i = 1:length(sublist)
     sub = sublist{i};
     for j = 1:length(sessions)
         sess = sessions{j};
         % set output directory
         mshbm_outdir = [mshbm_dir '/' sub '/' sess '/MSHBM_ts'];
         % find experience sampling .tsv files for fMRI
         es_tsv_files = [];
         es_tsv_files = dir([bids_dir '/' sub '/' sess '/func/sub*Exp*.tsv' ]);
         es_tsv_dir = es_tsv_files(1).folder;
         es_tsv_files = {es_tsv_files.name};
         % loop through ES runs
            for k = 1:length(es_tsv_files)
                data = []; header = []; raw = [];
                [data,header,raw] = tsvread([es_tsv_dir '/' es_tsv_files{k}]);
                % Get n TRs
                fmri = [];
                fmri = MRIread([es_tsv_dir '/' es_tsv_files{k}(1:end-10) 'bold.nii.gz']);
                TR_labels_blank = [];
                TR_labels_blank = NaN(fmri.nframes,1);
                TR = fmri.tr/1000;
                TR_labels_sec = []; starting_TR = 0;
                for l = 1:fmri.nframes  
                    TR_labels_sec = [TR_labels_sec; starting_TR];
                    starting_TR = starting_TR+TR;
                end
                % find TP segments in seconds
                tp_start_inds = []; tp_start_times = []; tp_end_inds = []; tp_end_times = [];
                tp_start_inds = contains(raw(:,5),'att');
                tp_end_inds = contains(raw(:,5),'conf');
                tp_start_times = data(tp_start_inds,1)+8; % add 8-second shift to account for HRF
                tp_end_times = data(tp_end_inds,3);
                % find rest segments in seconds
                rest_start_times = []; rest_end_times = [];
                rest_start_times = [tp_end_times]+8; % add 8-second shift after TP end to account for HRF
                rest_start_times(end)=[]; % remove last rest window (after last TP, because run ends soon)
                rest_start_times = [0; rest_start_times]; % add in first rest window start
                rest_end_times = tp_start_times-8; % unshift 8-seconds for rest offset
                % map TP segments to TRs and set TP labels
                TR_cond_labels = TR_labels_blank;
                for m = 1:length(tp_start_times)
                    [c,tp_start_TRs(m)] = min(abs(TR_labels_sec-tp_start_times(m)));
                    [c,tp_end_TRs(m)] = min(abs(TR_labels_sec-tp_end_times(m)));
                    TR_cond_labels(tp_start_TRs(m):tp_end_TRs(m)) = 2;
                end
                % map rest segments to TRs and set rest labels
                for m = 1:length(rest_start_times)
                    [c,rest_start_TRs(m)] = min(abs(TR_labels_sec-rest_start_times(m)));
                    [c,rest_end_TRs(m)] = min(abs(TR_labels_sec-rest_end_times(m)));
                    TR_cond_labels(rest_start_TRs(m):rest_end_TRs(m)) = 1;
                end
            % save TR condition labels for this run
            dlmwrite([mshbm_outdir '/ConditionLabel_' sub '_' sess '_' char(mshbm_es_run_labels(k)) '.txt'],TR_cond_labels);
            display(['done ' sub ' ' sess ' ' es_tsv_files{k} ' (' char(mshbm_es_run_labels(k)) ')']);
            end
     end
 end