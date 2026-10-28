% Extract Null DN time series from all sessions+runs+subjects
% settings
addpath(genpath('/Users/ak4379/Documents/MATLAB/cifti-matlab'));
% get list of subjects
sublist = dir('/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/CBIG_preproc_fmap_highpass/sub*');
sublist = {sublist.name};
sessions = {'ses-001'; 'ses-002'};
localizer = '300';

% defaults
CBIG_preproc_dir = '/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/CBIG_preproc_fmap_highpass'; % location of preprocessed fMRI data from CBIG 2016 preprocessing pipeline
%parpool(8);

% loops through all sessions and runs
for k=1:length(sublist)
for i = 1:length(sessions)
   % find all runs in this session
   lh_session_runs = dir([CBIG_preproc_dir '/' sublist{k} '/' sessions{i} '/' sublist{k} '/surf/lh*bld*']);
   %rh_session_runs = dir([CBIG_preproc_dir '/' sub '/' sessions{i} '/' sub '/surf/rh*bld*']);
   % extract time series for each run
   for j = 1:length(lh_session_runs)
       run_num = [];
       run_num = strfind(lh_session_runs(j).name,'bld');
       run_num = lh_session_runs(j).name(run_num+3:run_num+5);
       tic;
       ArealMSHBM_meants_highpass_NullDN(sublist{k},sessions{i},run_num,localizer);
       display(['done extracting areal MSHBM time series for ' sublist{k} ' ' sessions{i} ' run ' num2str(j)]);
       toc;
   end
end
end