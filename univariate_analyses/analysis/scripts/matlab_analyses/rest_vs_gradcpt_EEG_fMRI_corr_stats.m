% Compare EEG-fMRI correlations for each network: GradCPT vs. Rest
% Significant EEG-fMRI correlations for each network
% Must first run mk_EEG_fMRI_corr_matrices.m and mk_EEG_fMRI_corr_matrices_rest.m

beta_split = 1; % 0 = don't split up beta band; 1 = split up beta band
data_dir = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/correlation_data/corr_mats'];
mkdir('figs');
sublist = {'sub-001'; 'sub-002';'sub-004';'sub-006';'sub-007';'sub-008';'sub-009';'sub-010';...
    'sub-011';'sub-012';'sub-013';'sub-014';'sub-015';'sub-016'; 'sub-018';'sub-019';'sub-020';...
    'sub-021';'sub-023';'sub-024'};
sessions = {'ses-001';'ses-002'};
networks = {'DNa';'DNb';'FPCNa';'FPCNb';'dATNa';'dATNb'};
if beta_split==0
    freqs = {'δ'; 'θ'; 'α'; 'β'; 'γ'};
elseif beta_split==1
    freqs = {'δ'; 'θ'; 'α'; 'β1'; 'β2'; 'γ'};
end
skip_sess1 = {'sub-003';'sub-025'};
skip_sess2 = {'sub-017'};

colors = cbrewer('div','RdBu',256);
colors = abs(colors);
colors = flipud(colors);
inter_colors = cbrewer('div','PuOr',256);
inter_colors = abs(inter_colors);
inter_colors = flipud(inter_colors);

% Load in GradCPT matrices (runs 001) and average within subjects
for i=1:length(networks)
    for j=1:length(sublist)
        if ~(sum(contains(skip_sess2,sublist{j}))>0) && ~(sum(contains(skip_sess1,sublist{j}))>0)
        gradCPT_matrices = []; subject_mean = [];
        for k=1:length(sessions)
            corr_mat_canonical = [];
            if beta_split==0
                load([data_dir '/' sublist{j} '/' sessions{k} '/run-001_' networks{i} '_canonical.mat']);
            elseif beta_split==1
                load([data_dir '/' sublist{j} '/' sessions{k} '/run-001_' networks{i} '_canonical_beta.mat']);
                corr_mat_canonical = corr_mat_canonical_beta;
            end
            corr_mat_canonical(7:end,:) = []; % chop out lags > 10 sec
            corr_mat_canonical = fisherz(corr_mat_canonical); % fisher-z transform
            gradCPT_matrices(:,:,k) = corr_mat_canonical;
        end
        subject_mean = mean(gradCPT_matrices,3); 
        gradCPT_mean{i}(:,:,j) = subject_mean;
        else % for subjects with only one session
            if sum(contains(skip_sess1,sublist{j}))>0
                if beta_split==0
                    load([data_dir '/' sublist{j} '/ses-002/run-001_' networks{i} '_canonical.mat']);
                elseif beta_split==1
                    load([data_dir '/' sublist{j} '/ses-002/run-001_' networks{i} '_canonical_beta.mat']);
                    corr_mat_canonical = corr_mat_canonical_beta;
                end
            corr_mat_canonical(7:end,:) = []; % chop out lags > 10 sec
            corr_mat_canonical = fisherz(corr_mat_canonical); % fisher-z transform
            gradCPT_mean{i}(:,:,j) = corr_mat_canonical;  
            end

            if sum(contains(skip_sess2,sublist{j}))>0
                if beta_split==0
                    load([data_dir '/' sublist{j} '/ses-001/run-001_' networks{i} '_canonical.mat']);
                elseif beta_split==1
                    load([data_dir '/' sublist{j} '/ses-001/run-001_' networks{i} '_canonical_beta.mat']);
                    corr_mat_canonical = corr_mat_canonical_beta;
                end
            corr_mat_canonical(7:end,:) = []; % chop out lags > 10 sec
            corr_mat_canonical = fisherz(corr_mat_canonical); % fisher-z transform
            gradCPT_mean{i}(:,:,j) = corr_mat_canonical;  
            end

        end
    end
end

% Load in rest matrices (runs 002-004) and average within subjects
for i=1:length(networks)
    for j=1:length(sublist)
        if ~(sum(contains(skip_sess2,sublist{j}))>0) && ~(sum(contains(skip_sess1,sublist{j}))>0)
        rest_matrices = []; subject_mean = [];
        for k=1:length(sessions)
            rest_matrix =[]; rest_matrix_sess = []; rest_matrices =[];
            if beta_split==0
                rest_matrix = load([data_dir '/' sublist{j} '/' networks{i} '_mean_canonical_rest.mat']);
                rest_matrix.subject_mean_canonical(7:end,:) = []; % chop out lags > 10 sec
                rest_matrix = rest_matrix.subject_mean_canonical;
            elseif beta_split==1
                rest_matrix = load([data_dir '/' sublist{j} '/' networks{i} '_mean_canonical_beta_rest.mat']);
                rest_matrix.subject_mean_canonical_beta(7:end,:) = []; % chop out lags > 10 sec
                rest_matrix = rest_matrix.subject_mean_canonical_beta;
            end
            rest_matrix_sess(:,:,k) = rest_matrix;
        end
        subject_mean = mean(rest_matrix_sess,3); 
        rest_mean{i}(:,:,j) = subject_mean;

        else % for subjects with only one session
            if sum(contains(skip_sess1,sublist{j}))>0
                if beta_split==0
                    load([data_dir '/' sublist{j} '/ses-002/' networks{i} '_mean_canonical_rest.mat']);
                elseif beta_split==1
                    load([data_dir '/' sublist{j} '/ses-002/' networks{i} '_mean_canonical_beta_rest.mat']);
                    corr_mat_canonical = corr_mat_canonical_beta;
                end
            corr_mat_canonical(7:end,:) = []; % chop out lags > 10 sec
            corr_mat_canonical = fisherz(corr_mat_canonical); % fisher-z transform
            rest_mean{i}(:,:,j) = corr_mat_canonical;  
            end

            if sum(contains(skip_sess2,sublist{j}))>0
                if beta_split==0
                    load([data_dir '/' sublist{j} '/ses-001/' networks{i} '_mean_canonical_rest.mat']);
                elseif beta_split==1
                    load([data_dir '/' sublist{j} '/ses-001/' networks{i} '_mean_canonical_beta_rest.mat']);
                    corr_mat_canonical = corr_mat_canonical_beta;
                end
            corr_mat_canonical(7:end,:) = []; % chop out lags > 10 sec
            corr_mat_canonical = fisherz(corr_mat_canonical); % fisher-z transform
            rest_mean{i}(:,:,j) = corr_mat_canonical;  
            end
        end
    end
end

% similarity between GradCPT and Rest heatmaps per subject
for i=1:length(networks)
    for j=1:length(sublist)
        curr_rest = []; curr_gradcpt = [];
        curr_rest = rest_mean{i}(:,:,j);
        curr_gradcpt = gradCPT_mean{i}(:,:,j);
        r_intertask(i,j) = corr(curr_rest(:),curr_gradcpt(:),'Type','Pearson');
    end
end

% Loop through networks and perform stats for rest and GradCPT
for i=1:length(networks)
    rest_network_allsubs = []; rest_p_mask = []; rest_subject_mean = [];
    gradcpt_network_allsubs = []; gradcpt_p_mask = []; gradcpt_subject_mean = [];
    rest_network_allsubs = rest_mean{i};
    gradcpt_network_allsubs = gradCPT_mean{i};
    rest_allnetworks_allsubs{i} = rest_network_allsubs; % store all networks
    gradcpt_allnetworks_allsubs{i} = gradcpt_network_allsubs;
    % mass univariate stats for each network
    rest_p_mat = NaN(size(rest_network_allsubs,1),size(rest_network_allsubs,2));
    gradcpt_p_mat = NaN(size(gradcpt_network_allsubs,1),size(gradcpt_network_allsubs,2));
    for j = 1:size(rest_network_allsubs,1) % for each lag
       for k = 1:size(rest_network_allsubs,2) % for each frequency band
           rest_p_mat(j,k) = signtest(squeeze(rest_network_allsubs(j,k,:)));
           gradcpt_p_mat(j,k) = signtest(squeeze(gradcpt_network_allsubs(j,k,:)));
       end
    end
    % fdr correction
    [rest_p_fdr,rest_p_mask] = fdr(rest_p_mat,.05);
    [gradcpt_p_fdr,gradcpt_p_mask] = fdr(gradcpt_p_mat,.05);
    % store mean and p values across subjects for plotting later
    rest_network_mean{i} = mean(rest_network_allsubs,3);
    gradcpt_network_mean{i} = mean(gradcpt_network_allsubs,3);
    rest_network_p{i} = rest_p_mask;
    gradcpt_network_p{i} = gradcpt_p_mask;
end

% paired tests between GradCPT vs Rest for each network
for i=1:length(networks)
p_fdr=[]; intertask_p_mask =[];
    %network_pair = network_pairs{i};
    %network_a = allnetworks_allsubs{network_pair(1)};
    %network_b = allnetworks_allsubs{network_pair(2)};
    %signrank
    for j=1:size(rest_mean{i},1)
        for k=1:size(rest_mean{i},2)
            [intertask_p{i}(j,k),h,z] = signrank(squeeze(gradCPT_mean{i}(j,k,:)),squeeze(rest_mean{i}(j,k,:)));
            intertask_z{i}(j,k) = z.zval;
        end
    end
    % fdr correction
    [p_fdr,intertask_p_mask] = fdr(intertask_p{i},.05);
    intertask_p_fdr{i} = intertask_p_mask;
end

% Rest: plot mean heatmaps across all subjects for each network
figure('Position',[200,200,215,1800]);
for i=1:length(rest_network_mean)
    subplot(6,1,i)
    imagesc(rest_network_mean{i}); % Display the image
    colormap(colors);
    set(gca, 'YDir', 'normal') % flip y 
    c=colorbar('Location','eastoutside');
    c.FontSize=10;
    c.Label.String=[{['z (\rho)' '_{EEG,fMRI}']}];
    %c.Label.Position = [4 0];
    set(gca,'xtick',1:size(rest_network_mean{i},2))
    set(gca,'ytick',1:size(rest_network_mean{i},1))
    set(gca,'xticklabel',freqs);
    set(gca,'ytick',1:size(rest_network_mean{i},1))
    set(gca,'yticklabel',0:2:10)
    set(gca,'TickLength',[0 0])
    set(gca,'FontSize',12,'FontName','Arial');
    %xlabel({['Frequency band'];[' ']});
    ylabel({['Lag relative '];['to BOLD (s)']});
    c_max = max(abs(rest_network_mean{i}(:)));
    %caxis([-.04 .04]);
    caxis([-c_max c_max]);
    title([networks{i}],'Fontweight','normal')
    box off
    hold on;
    [y_fdr,x_fdr]=find(rest_network_p{i}==1);
    s=scatter(x_fdr,y_fdr,60,'k','*','MarkerEdgeColor',[255 215 0]/255);
    hold on;
end
print('-opengl','-r600','-dpng',['figs/rest_heatmaps.png']);  
pause; close;

% GradCPT: plot mean heatmaps across all subjects for each network
figure('Position',[200,200,215,1800]);
for i=1:length(gradcpt_network_mean)
    subplot(6,1,i)
    imagesc(gradcpt_network_mean{i}); % Display the image
    colormap(colors);
    set(gca, 'YDir', 'normal') % flip y 
    c=colorbar('Location','eastoutside');
    c.FontSize=10;
    c.Label.String=[{['z (\rho)' '_{EEG,fMRI}']}];
    %c.Label.Position = [4 0];
    set(gca,'xtick',1:size(gradcpt_network_mean{i},2))
    set(gca,'ytick',1:size(gradcpt_network_mean{i},1))
    set(gca,'xticklabel',freqs);
    set(gca,'ytick',1:size(gradcpt_network_mean{i},1))
    set(gca,'yticklabel',0:2:10)
    set(gca,'TickLength',[0 0])
    set(gca,'FontSize',12,'FontName','Arial');
    %xlabel({['Frequency band'];[' ']});
    ylabel({['Lag relative '];['to BOLD (s)']});
    c_max = max(abs(gradcpt_network_mean{i}(:)));
    caxis([-c_max c_max]);
    %caxis([-.04 .04]);
    title([networks{i}],'Fontweight','normal')
    box off
    hold on;
    [y_fdr,x_fdr]=find(gradcpt_network_p{i}==1);
    s=scatter(x_fdr,y_fdr,60,'k','*','MarkerEdgeColor',[255 215 0]/255);
    hold on;
end
print('-opengl','-r600','-dpng',['figs/gradcpt_heatmaps.png']);  
pause; close;

% plot GradCPT vs. Rest
figure('Position',[200,200,215,1800]);
for i=1:length(intertask_z)
    subplot(6,1,i)
    imagesc(intertask_z{i}); % Display the image
    colormap(inter_colors);
    set(gca, 'YDir', 'normal') % flip y 
    c=colorbar('Location','eastoutside');
    c.FontSize=14;
    c.Label.String=[{['z-score']}];
    %c.Label.Position = [4 0];
    set(gca,'xtick',1:size(intertask_z{i},2))
    set(gca,'ytick',1:size(intertask_z{i},1))
    set(gca,'xticklabel',freqs);
    set(gca,'ytick',1:size(intertask_z{i},1))
    set(gca,'yticklabel',0:2:10)
    set(gca,'TickLength',[0 0])
    set(gca,'FontSize',12,'FontName','Arial');
    %xlabel({['Frequency band'];[' ']});
    ylabel({['Lag relative '];['to BOLD (s)']});
    c_max = max(abs(intertask_z{i}(:)));
    caxis([-c_max c_max]);
    title([networks{i}],'Fontweight','normal')
    box off
    hold on;
    [y_fdr,x_fdr]=find(intertask_p_fdr{i}==1);
    s=scatter(x_fdr,y_fdr,60,'k','*','MarkerEdgeColor',[255 215 0]/255);
    hold on;
end
print('-opengl','-r600','-dpng',['figs/intertask_heatmaps.png']); 
pause; close;


