% Run mass-univariate stats and plot
% Significant EEG-fMRI correlations for each network
% Significant differences in EEG-fMRI correlations between pairs of subnetworks
% Must first run 01_mk_EEG_fMRI_corr_matrices.m

data_dir = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/correlation_data/corr_mats'];
mkdir('figs');
sublist = {'sub-001'; 'sub-002';'sub-004';'sub-006';'sub-007';'sub-008';'sub-009';'sub-010';...
    'sub-011';'sub-012';'sub-013';'sub-014';'sub-015';'sub-016';'sub-018';'sub-019';'sub-020';...
    'sub-021';'sub-023';'sub-024'};
sessions = {'ses-001';'ses-002'};
networks = {'DNa';'DNb';'FPCNa';'FPCNb';'dATNa';'dATNb'};
network_pairs = {[1 2];[3 4];[5 6]};
pair_names = {'DNa > DNb'; 'FPCNa > FPCNb'; 'dATNa > dATNb'};
freqs = {'δ'; 'θ'; 'α'; 'β'; 'γ'};

colors = cbrewer('div','RdBu',256);
colors = abs(colors);
colors = flipud(colors);
inter_colors = cbrewer('div','PuOr',256);
inter_colors = abs(inter_colors);
inter_colors = flipud(inter_colors);

% Load in gradCPT matrices and average within subjects
for i=1:length(networks)
    for j=1:length(sublist)
        gradCPT_matrices = []; subject_mean = [];
        for k=1:length(sessions)
            corr_mat_canonical = [];
            load([data_dir '/' sublist{j} '/' sessions{k} '/run-001_' networks{i} '_canonical.mat']);
            corr_mat_canonical(7:end,:) = []; % chop out lags > 10 sec
            corr_mat_canonical = fisherz(corr_mat_canonical); % fisher-z transform
            gradCPT_matrices(:,:,k) = corr_mat_canonical;
        end
        subject_mean = mean(gradCPT_matrices,3); 
        gradCPT_mean{i}(:,:,j) = subject_mean;
    end
end

% Loop through networks and perform stats
for i=1:length(networks)
    network_allsubs = []; p_mask = [];
    subject_mean = [];
    network_allsubs = gradCPT_mean{i};
    allnetworks_allsubs{i} = network_allsubs; % store all networks
    % mass univariate stats for each network
    p_mat = NaN(size(network_allsubs,1),size(network_allsubs,2));
    for j = 1:size(network_allsubs,1) % for each lag
       for k = 1:size(network_allsubs,2) % for each frequency band
           p_mat(j,k) = signtest(squeeze(network_allsubs(j,k,:)));
       end
    end
    % fdr correction
    [p_fdr,p_mask] = fdr(p_mat,.05);
    % store mean and p values across subjects for plotting later
    network_mean{i} = mean(network_allsubs,3);
    network_p{i} = p_mask;
end

% paired tests between subnetworks
for i = 1:size(allnetworks_allsubs,2)/2
    p_fdr=[]; p_mask =[];
    network_pair = network_pairs{i};
    network_a = allnetworks_allsubs{network_pair(1)};
    network_b = allnetworks_allsubs{network_pair(2)};
    %signrank
    for j=1:size(network_a,1)
        for k=1:size(network_a,2)
            [internetwork_p{i}(j,k),h,z] = signrank(squeeze(network_a(j,k,:)),squeeze(network_b(j,k,:)));
            internetwork_z{i}(j,k) = z.zval;
        end
    end
    % fdr correction
    [p_fdr,p_mask] = fdr(internetwork_p{i},.05);
    internetwork_p_fdr{i} = p_mask;
end

% plot mean heatmaps across all subjects for each network
figure('Position',[200,200,1500,150]);
for i=1:length(network_mean)
    subplot(1,6,i)
    imagesc(network_mean{i}); % Display the image
    colormap(colors);
    set(gca, 'YDir', 'normal') % flip y 
    c=colorbar('Location','eastoutside');
    c.FontSize=10;
    c.Label.String=['z (\rho)' '_{EEG,fMRI}'];
    set(gca,'xtick',1:size(network_mean{i},2))
    set(gca,'ytick',1:size(network_mean{i},1))
    set(gca,'xticklabel',freqs);
    set(gca,'ytick',1:size(network_mean{i},1))
    set(gca,'yticklabel',0:2:10)
    set(gca,'TickLength',[0 0])
    set(gca,'FontSize',12,'FontName','Arial');
    xlabel({['Frequency band'];[' ']});
    if i~=1
        ylabel({[' '];[' ']});
    else
        ylabel({['Lag relative '];['to BOLD (s)']});
    end
    c_max = max(abs(network_mean{i}(:)));
    caxis([-c_max c_max]);
    title(networks{i})
    box off
    hold on;
    [y_fdr,x_fdr]=find(network_p{i}==1);
    s=scatter(x_fdr,y_fdr,60,'k','*','MarkerEdgeColor',[255 215 0]/255);
    hold on;
end
print('-opengl','-r600','-dpng',['figs/gradCPT_heatmaps.png']);  
pause; close;

% plot network a vs b difference maps
figure('Position',[200,200,300,900]);
for i=1:length(internetwork_z)
    subplot(3,1,i)
    imagesc(internetwork_z{i}); % Display the image
    colormap(inter_colors);
    set(gca, 'YDir', 'normal') % flip y 
    c=colorbar('Location','eastoutside');
    c.FontSize=14;
    c.Label.String=['z-score'];
    set(gca,'xtick',1:size(internetwork_z{i},2))
    set(gca,'ytick',1:size(internetwork_z{i},1))
    set(gca,'xticklabel',freqs);
    set(gca,'ytick',1:size(internetwork_z{i},1))
    set(gca,'yticklabel',0:2:10)
    set(gca,'TickLength',[0 0])
    set(gca,'FontSize',16,'FontName','Arial');
    xlabel({['Frequency band'];[' ']});
    ylabel({[' '];[' '];['Lag relative to BOLD (s)']});
    c_max = max(abs(internetwork_z{i}(:)));
    caxis([-c_max c_max]);
    title([pair_names{i}])
    box off
    hold on;
    [y_fdr,x_fdr]=find(internetwork_p_fdr{i}==1);
    s=scatter(x_fdr,y_fdr,60,'k','*','MarkerEdgeColor',[255 215 0]/255);
    hold on;
end
print('-opengl','-r600','-dpng',['figs/gradCPT_internetwork_heatmaps.png']); 
pause; close;


