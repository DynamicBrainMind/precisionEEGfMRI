function extract_NullDN_BOLD_ts(sub,sess,run_num,localizer)

% extract time series from DN permutations (generated from step1_mk_DN_null_labels.m)
% Inputs: 
% 1. subject
% 2. session
% 3. run number (as string; e.g. '001')
% 4. name of MSHBM localizer folder (e.g. 'Localizer300'), which should indicate n parcels

% Defaults
n_perm = 1000; % number of surrogate DN's to generate 
CBIG_preproc_dir = '/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/CBIG_preproc_fmap_highpass'; % location of data preprocessed with CBIG 2016 preprocessing pipeline
MSHBM_dir = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/ArealMSHBM_sdc/' sub '/' sess '/' localizer '/ind_parcellation_gMSHBM/test_set/4_sess/beta5/NullDN']; 
outdir = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/ArealMSHBM_sdc/' sub '/' sess '/' localizer '/ArealMSHBM_ts/NullDN'];
[~,~]=mkdir([outdir]);

% Load in fmri data in fsaverage6 space (from CBIG pipeline)
surf_files = dir([CBIG_preproc_dir '/' sub '/' sess '/' sub '/surf/*bld' run_num '*']);
lh_file = surf_files(1).name;
rh_file = surf_files(2).name;
lh_fmri = MRIread([CBIG_preproc_dir '/' sub '/' sess '/' sub '/surf/' surf_files(1).name]);
lh_vol = lh_fmri.vol;
lh_vol_size = size(lh_vol);
lh_vol = reshape(lh_vol, prod(lh_vol_size(1:3)), lh_vol_size(4));
rh_fmri = MRIread([CBIG_preproc_dir '/' sub '/' sess '/' sub '/surf/' surf_files(2).name]);
rh_vol = rh_fmri.vol;
rh_vol_size = size(rh_vol);
rh_vol = reshape(rh_vol, prod(rh_vol_size(1:3)), rh_vol_size(4));

% Load motion outliers (censor list)
run_num = lh_file(strfind(lh_file,'bld'):strfind(lh_file,'bld')+5);
mc_file = dir([CBIG_preproc_dir '/' sub '/' sess '/' sub '/qc/*' run_num '_FDRMS*motion_outliers.txt']);
censor_list=load([mc_file.folder '/' mc_file.name]);
display(['censor file is ' mc_file.name]);

%% Loop through null DN permutations
% Load in MSHBM-derived areal labels
%parfor p = 1:n_perm 
for p = 1:n_perm 
    %display(['Permutation number: ' num2str(p)]);
    data = []; lh_labels = []; rh_labels = [];
    data = cifti_read ([MSHBM_dir '/' sub '_NullDNs_' num2str(p) '.dlabel.nii']);
    lh_labels = data.cdata(1:40962);
    rh_labels = data.cdata(40963:end);

    % extract time series from each network
    unique_lh = unique(lh_labels); 
    unique_lh(1) = []; 
    unique_rh = unique(rh_labels); 
    unique_rh(1) = []; 

    % extract time series from all MSHBM-derived networks
    allnetwork_ts = [];
    for i = 1:max(lh_labels)
        lh = []; rh = []; network_ts = [];
        lh = lh_vol(lh_labels==i,:);
        rh = rh_vol(rh_labels==i,:);
        network_ts = [lh; rh];
        allnetwork_ts = [allnetwork_ts mean(network_ts,1)'];
    end
    % get NullDNa and NullDNb
    allnetwork_ts = allnetwork_ts(:,301:302);

    % set censored time points to NaN
    allnetwork_ts(censor_list==0,:)=NaN;

    % save time series as text file
    dlmwrite([outdir filesep sess '_' run_num '_NullDNa_highpass_' num2str(p) '.txt'],allnetwork_ts(:,1));
    dlmwrite([outdir filesep sess '_' run_num '_NullDNb_highpass_' num2str(p) '.txt'],allnetwork_ts(:,2));
end

