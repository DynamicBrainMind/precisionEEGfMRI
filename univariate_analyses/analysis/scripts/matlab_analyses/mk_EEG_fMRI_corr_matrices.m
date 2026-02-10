% set output location
outdir = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/correlation_data'];
mkdir([outdir '/corr_mats']);
skip_sess1 = {'sub-003';'sub-025'};
skip_sess2 = {'sub-017'};

% Load correlations (averaged across channels)
data = readtable([outdir '/correlations_long_across_channel.csv']);
%textdata = data.textdata;
corr_data = data.cors;
%corr_data = [NaN; corr_data]; % shift by 1 to match textdata
subs = unique(data.subject);
%subs = unique(data.textdata(2:end,1));
sessions = unique(data.session);
%sessions = unique(data.textdata(2:end,2));
runs = unique(data.run);
%runs = unique(data.textdata(2:end,3));
freqs = unique(data.frequency);
%freqs = unique(data.textdata(2:end,4));
lags = unique(data.lag);
%lags = unique(data.textdata(2:end,5));
networks = unique(data.network);
%networks = unique(data.textdata(2:end,6));
textdata = table2cell(data);

% Compile mean corr matrices for each run (one for each network)
lags = sort(lags);
freqs = sort(freqs);

for i = 1:length(subs)
    mkdir([outdir '/corr_mats/' subs{i}]);
    for j = 1:length(sessions)
        corr_mat_allnetworks = []; corr_mat_canonical_allnetworks = []; corr_mat_canonical_beta_allnetworks = [];
    if ~(sum(contains(skip_sess2,subs{i}))>0 && j==2) && ~(sum(contains(skip_sess1,subs{i}))>0 && j==1)
    mkdir([outdir '/corr_mats/' subs{i} '/' sessions{j}]);
        for k = 1:length(runs)
            for m = 1:length(networks)
                curr_textdata = []; curr_data = []; sub_sess_run_network_inds = []; sub_inds = [];
                sess_inds = []; run_inds = []; network_inds = []; sub_sess_inds = []; sub_sess_run_inds = [];
                corr_mat_canonical = []; corr_mat_canonical_beta = [];
                sub_inds = strmatch(subs{i},data.subject,'exact');
                %sub_inds = strmatch(subs{i},textdata(:,1),'exact');
                sess_inds = strmatch(sessions{j},data.session,'exact');
                %sess_inds = strmatch(sessions{j},textdata(:,2),'exact'); 
                run_inds = strmatch(runs{k},data.run,'exact');
                %run_inds = strmatch(runs{k},textdata(:,3),'exact');
                network_inds = strmatch(networks{m},data.network,'exact');
                %network_inds = strmatch(networks{m},textdata(:,6),'exact');
                sub_sess_inds = intersect(sub_inds,sess_inds);
                sub_sess_run_inds = intersect(sub_sess_inds,run_inds);
                sub_sess_run_network_inds = intersect(sub_sess_run_inds,network_inds);
                curr_textdata = textdata(sub_sess_run_network_inds,:);
                curr_data = corr_data(sub_sess_run_network_inds);
                % fill in matrix (x = freq, y = lag)
                corr_mat = NaN(length(lags),length(freqs));
                curr_lags = cell2mat(curr_textdata(:,5));
                curr_freqs = cell2mat(curr_textdata(:,4));
                for n = 1:length(lags)
                    for p = 1:length(freqs)
                        curr_lag = lags(n);
                        lag_inds = find(curr_lags==curr_lag);
                        freq_inds = find(curr_freqs==p);
                        curr_ind = intersect(lag_inds,freq_inds);
                        corr_mat(n,p) = curr_data(curr_ind);                      
                    end 
                end
                % make corr matrix with averaging of canonical frequency bands
                corr_mat_canonical = [mean(corr_mat(:,1:3),2) mean(corr_mat(:,4:7),2) mean(corr_mat(:,8:12),2)...
                    mean(corr_mat(:,13:30),2) mean(corr_mat(:,31:40),2)];
                % make corr matrix with averaging of canonical frequency: with beta subbands included
                %corr_mat_canonical_beta = [mean(corr_mat(:,1:3),2) mean(corr_mat(:,4:7),2) mean(corr_mat(:,8:12),2)...
                %    mean(corr_mat(:,13:16),2),mean(corr_mat(:,17:23),2),mean(corr_mat(:,24:30),2) , mean(corr_mat(:,31:40),2)];
                corr_mat_canonical_beta = [mean(corr_mat(:,1:3),2) mean(corr_mat(:,4:7),2) mean(corr_mat(:,8:12),2)...
                    mean(corr_mat(:,13:20),2),mean(corr_mat(:,21:30),2), mean(corr_mat(:,31:40),2)];
                % store all correlations
                corr_mat_allnetworks{k}(:,:,m) = corr_mat;
                corr_mat_canonical_allnetworks{k}(:,:,m) = corr_mat_canonical;
                corr_mat_canonical_beta_allnetworks{k}(:,:,m) = corr_mat_canonical_beta;
                save([outdir '/corr_mats/' subs{i} '/' sessions{j} '/' runs{k} '_' networks{m} '.mat'],'corr_mat');
                save([outdir '/corr_mats/' subs{i} '/' sessions{j} '/' runs{k} '_' networks{m} '_canonical.mat'],'corr_mat_canonical');
                save([outdir '/corr_mats/' subs{i} '/' sessions{j} '/' runs{k} '_' networks{m} '_canonical_beta.mat'],'corr_mat_canonical_beta');
                display(['done ' subs{i} ' ' sessions{j} ' ' runs{k} ' ' networks{m}]);
            end
        end
    
    %Compute and save average corr matrices for sessions 
    for k = 1:length(networks)
        network_allruns = []; session_mean = []; network_allruns_canonical = []; network_allruns_canonical_beta = []; 
        for m = 1:size(corr_mat_allnetworks,2)
            network_allruns(:,:,m) = corr_mat_allnetworks{m}(:,:,k);
            network_allruns_canonical(:,:,m) = corr_mat_canonical_allnetworks{m}(:,:,k);
            network_allruns_canonical_beta(:,:,m) = corr_mat_canonical_beta_allnetworks{m}(:,:,k);
        end
        session_mean = mean(network_allruns,3);
        session_mean_canonical = mean(network_allruns_canonical,3);
        session_mean_canonical_beta = mean(network_allruns_canonical_beta,3);
        save([outdir '/corr_mats/' subs{i} '/' sessions{j} '/' networks{k} '_mean.mat'],'session_mean');
        save([outdir '/corr_mats/' subs{i} '/' sessions{j} '/' networks{k} '_mean_canonical.mat'],'session_mean_canonical');
        save([outdir '/corr_mats/' subs{i} '/' sessions{j} '/' networks{k} '_mean_canonical_beta.mat'],'session_mean_canonical_beta');
    end
    end
    end
end

 % Compute average corr matrices for subject
 for i=1:length(subs)
     if isempty(strmatch(subs{i},[skip_sess1; skip_sess2],'exact'))
     for j=1:length(networks)
         network_mat = []; subject_mean = []; network_mat_canonical = []; network_mat_canonical_beta = [];
         for k=1:length(sessions)
           network_corr = load([outdir '/corr_mats/' subs{i} '/' sessions{k} '/' networks{j} '_mean.mat']);
           network_mat(:,:,k) = network_corr.session_mean;
           network_corr_canonical = load([outdir '/corr_mats/' subs{i} '/' sessions{k} '/' networks{j} '_mean_canonical.mat']);
           network_mat_canonical(:,:,k) = network_corr_canonical.session_mean_canonical;
           network_corr_canonical_beta = load([outdir '/corr_mats/' subs{i} '/' sessions{k} '/' networks{j} '_mean_canonical_beta.mat']);
           network_mat_canonical_beta(:,:,k) = network_corr_canonical_beta.session_mean_canonical_beta;
         end
         subject_mean = mean(network_mat,3);
         subject_mean_canonical = mean(network_mat_canonical,3);
         subject_mean_canonical_beta = mean(network_mat_canonical_beta,3);
         save([outdir '/corr_mats/' subs{i} '/' networks{j} '_mean.mat'],'subject_mean');
         save([outdir '/corr_mats/' subs{i} '/' networks{j} '_mean_canonical.mat'],'subject_mean_canonical');
         save([outdir '/corr_mats/' subs{i} '/' networks{j} '_mean_canonical_beta.mat'],'subject_mean_canonical_beta');
     end
     end 
 end
 
 % for subjects with only 1 session, copy session mean to be subject mean
 for i=1:length(skip_sess1)
     for j=1:length(networks)
        load([[outdir '/corr_mats/' skip_sess1{i} '/ses-002/' networks{j} '_mean.mat']]);
        subject_mean = session_mean;
        load([[outdir '/corr_mats/' skip_sess1{i} '/ses-002/' networks{j} '_mean_canonical.mat']]);
        subject_mean_canonical = session_mean_canonical;    
        load([[outdir '/corr_mats/' skip_sess1{i} '/ses-002/' networks{j} '_mean_canonical_beta.mat']]);
        subject_mean_canonical_beta = session_mean_canonical_beta;
        save([outdir '/corr_mats/' skip_sess1{i} '/' networks{j} '_mean.mat'],'subject_mean');
        save([outdir '/corr_mats/' skip_sess1{i} '/' networks{j} '_mean_canonical.mat'],'subject_mean_canonical');
        save([outdir '/corr_mats/' skip_sess1{i} '/' networks{j} '_mean_canonical_beta.mat'],'subject_mean_canonical_beta');
     end
 end
 
  for i=1:length(skip_sess2)
     for j=1:length(networks)
        load([[outdir '/corr_mats/' skip_sess2{i} '/ses-001/' networks{j} '_mean.mat']]);
        subject_mean = session_mean;
        load([[outdir '/corr_mats/' skip_sess2{i} '/ses-001/' networks{j} '_mean_canonical.mat']]);
        subject_mean_canonical = session_mean_canonical;
        load([[outdir '/corr_mats/' skip_sess2{i} '/ses-001/' networks{j} '_mean_canonical_beta.mat']]);
        subject_mean_canonical_beta = session_mean_canonical_beta;
        save([outdir '/corr_mats/' skip_sess2{i} '/' networks{j} '_mean.mat'],'subject_mean');
        save([outdir '/corr_mats/' skip_sess2{i} '/' networks{j} '_mean_canonical.mat'],'subject_mean_canonical');
        save([outdir '/corr_mats/' skip_sess2{i} '/' networks{j} '_mean_canonical_beta.mat'],'subject_mean_canonical_beta');
     end
 end
 
