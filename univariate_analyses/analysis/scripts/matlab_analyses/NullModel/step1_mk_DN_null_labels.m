% Take DNa and DNb regions from areal MSHBM, and mix labels of regions to create a null DNa and DNb

%% get list of subjects
sublist = dir('/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/CBIG_preproc_fmap_highpass/sub*'); % get list from a directory containing subject folders
sublist = {sublist.name};
sessions = {'ses-001'; 'ses-002'};
nparcels = '300';
n_perm = 1000; % number of surrogate DN's to generate 

%% set paths
MSHBM_outdir='/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/ArealMSHBM_sdc'; % output directory location
addpath(genpath('/Users/ak4379/Documents/MATLAB/cifti-matlab')); % cifti-matlab repo location
addpath(genpath('/Users/ak4379/Documents/scripts/CBIG')); % CBIG repo location
n_sess = '4_sess'; % Name of subdirectory under "test_set" within each subject's ArealMSHBM output

%% Load a sample subject to set up parcel label randomization
data_template = cifti_read ([MSHBM_outdir '/sub-001/ses-001/' nparcels '/ind_parcellation_gMSHBM/test_set/' n_sess '/beta5/sub-001_fsaverage6_MSHBM_parcellation.dlabel.nii']);
% SET THE LINE ABOVE TO LOCATION OF PARCELLATION OUTPUTS

%% Loop through permuations
ari_all = NaN(n_perm,1);
parpool(8);
parfor p = 1:n_perm
%for p = 1:n_perm
    display(['Permutation number: ' num2str(p)]);
    % reset variables
    DNa_parcels=[]; DNb_parcels=[]; DNa_parcels_resample=[]; DNb_parcels_resample=[];
    % get DNa and DNb parcel numbers
    DNa_parcels = find(contains({data_template.diminfo{1,2}.maps.table.name},'DefaultA')==1);
    DNb_parcels = find(contains({data_template.diminfo{1,2}.maps.table.name},'DefaultB')==1);
    % concatenate and scramble these parcel numbers
    DN_parcels = [DNa_parcels DNb_parcels];
    orig_DN_labels = [ones(length(DNa_parcels),1); ones(length(DNb_parcels),1)*2];
    idx = randperm(length(DN_parcels));
    shuffled_DN_parcels = DN_parcels(idx);
    shuffled_DN_labels = orig_DN_labels(idx);
    % compute ARI between this permutation and original
    [ri,ari] = randindex(orig_DN_labels,shuffled_DN_labels);
    ari_all(p) = ari;
    % Label a random half of DNa parcels as DNb and vice versa, keeping
    % total number of parcels consistent with original data
    orig_labels = data_template.cdata;
    DNa_parcels_resample = shuffled_DN_parcels(1:length(DNa_parcels));
    DNb_parcels_resample = shuffled_DN_parcels(length(DNa_parcels)+1:end);

    %% Loop through subjects and sessions
    for i=1:length(sublist)
        for j=1:length(sessions)
            data = [];
            [~,~]=mkdir([MSHBM_outdir '/' sublist{i} '/' sessions{j} '/' nparcels '/ind_parcellation_gMSHBM/test_set/' n_sess '/beta5/NullDN'])
            output_dir=[MSHBM_outdir filesep sublist{i} filesep sessions{j} filesep nparcels];
            data = cifti_read ([MSHBM_outdir '/' sublist{i} '/' sessions{j} '/' nparcels '/ind_parcellation_gMSHBM/test_set/' n_sess...
                '/beta5/' sublist{i} '_fsaverage6_MSHBM_parcellation.dlabel.nii']);

        % Add new labels for null DNa and DNb that combine all parcels (including relabeled ones)
            data.diminfo{1,2}.maps.table(302).name = 'NullDNa';
            data.diminfo{1,2}.maps.table(302).key = 301;
            data.diminfo{1,2}.maps.table(302).rgba = data.diminfo{1,2}.maps.table(2).rgba;
            data.diminfo{1,2}.maps.table(303).name = 'NullDNb';
            data.diminfo{1,2}.maps.table(303).key = 302;
            data.diminfo{1,2}.maps.table(303).rgba = data.diminfo{1,2}.maps.table(13).rgba;
            for k = 1:length(DNa_parcels)
                data.cdata(find(data.cdata==DNa_parcels_resample(k)))=301;
            end
            for k = 1:length(DNb_parcels)
                data.cdata(find(data.cdata==DNb_parcels_resample(k)))=302;
            end

        % Remove all other labels (including removed DNa labels)
            data.cdata(find(data.cdata~=301 & data.cdata~=302))=0;

        % write out the permuted labels
            cifti_write(data,[MSHBM_outdir '/' sublist{i} '/' sessions{j} '/' nparcels '/ind_parcellation_gMSHBM/test_set/' n_sess '/beta5/NullDN/' sublist{i} '_NullDNs_' num2str(p) '.dlabel.nii']);
            display(['done ' sublist{i} ' ' sessions{j}]);
            
         % save the original merged DNa and DNb from areal MSHBM with other labels removed
        end
    end
end

% save ARI values
mkdir('results')
save(['results/ari.mat'],'ari_all');


    