% plot EEG-fMRI corrs for all frequency and up to 20 seconds of lag
% Must first run 01_mk_EEG_fMRI_corr_matrices.m

data_dir = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/correlation_data/corr_mats'];
mkdir('figs');
sublist = dir([data_dir '/sub*']);
sublist = {sublist.name};
%networks = {'Yeo7_DN'; 'DNa'; 'DNb'};
%networks = {'DNa';'DNb';'FPCNa';'FPCNb';'dATNa';'dATNb'};
%networks = {'RH_DefaultA_IPL'; 'RH_DefaultB_IPL'; 'RH_ContB_IPL'};
%networks = {'LH_DefaultA_PFCm'; 'LH_DefaultB_PFCd'; 'LH_ContB_PFCmp'};
%networks = {'LH_DefaultA_IPL'; 'LH_DefaultB_IPL'; 'LH_ContB_IPL'};
networks = {'RH_DefaultA_Temp'; 'RH_DefaultB_Temp'; 'RH_ContB_Temp'};
%networks = {'RH_DefaultA_pCunPCC'; 'RH_ContB_PCC'};
freqs = {'δ'; 'θ'; 'α'; 'β'; 'γ'};
do_network_pairs = 0; % set to 1 to do stats/plot network pairs
network_pairs = {[1 2];[3 4];[5 6]};
pair_names = {'DNa > DNb'; 'FPCNa > FPCNb'; 'dATNa > dATNb'};

colors = cbrewer('div','RdBu',256);
colors = abs(colors);
colors = flipud(colors);
inter_colors = cbrewer('div','PuOr',256);
inter_colors = abs(inter_colors);
inter_colors = flipud(inter_colors);

% Load in and concatenate matrices for each network
for i=1:length(networks)
    network_allsubs = []; p_mask = [];
    for j=1:length(sublist)
        subject_mean = [];
        load([data_dir '/' sublist{j} '/' networks{i} '_mean.mat']);
        %subject_mean_canonical(7:end,:) = []; % chop out lags > 10 sec
        %subject_mean = fisherz(subject_mean);
        network_allsubs(:,:,j) = subject_mean;
        allnetworks_allsubs{i} = network_allsubs; % store all networks
    end
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
if do_network_pairs==1
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
end

% plot mean heatmaps across all subjects for each network
figure('Position',[200,200,400,1500]);
for i=1:length(network_mean)
    subplot(3,1,i)
    imagesc(network_mean{i}); % Display the image
    colormap(colors);
    set(gca, 'YDir', 'normal') % flip y 
    c=colorbar('Location','eastoutside');
    c.FontSize=16;
    c.Label.String=[{['\rho' '_{EEG,fMRI}'];[' ']}];
    c.Label.VerticalAlignment = 'bottom';
    set(gca,'xtick',0:5:size(network_mean{i},2))
    set(gca,'ytick',1:size(network_mean{i},1))
    %set(gca,'xticklabel',1:5:40);
    set(gca,'ytick',1:size(network_mean{i},1))
    set(gca,'yticklabel',0:2:20)
    set(gca,'TickLength',[0 0])
    set(gca,'FontSize',18,'FontName','Arial');
    xlabel({['Frequency'];[' ']});
    %if mod(i, 2) == 0
    %    ylabel({[' '];[' '];[' ']});
    %else
        ylabel(['Lag relative to BOLD (s)']);
    %end
    c_max = max(abs(network_mean{i}(:)));
    caxis([-c_max c_max]);
        c_limits = c.Limits;
    c.Ticks = c_limits;
    %title(networks{i})
    box off
    hold on;
    [y_fdr,x_fdr]=find(network_p{i}==1);
    %s=scatter(x_fdr,y_fdr,60,'k','*','MarkerEdgeColor',[255 215 0]/255);
    hold on;
end
print('-opengl','-r600','-dpng',['figs/long_heatmaps.png']);  
pause; close;


