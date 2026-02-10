% Run mass-univariate stats and plot
% Significant EEG-fMRI correlations for each network
% Significant differences in EEG-fMRI correlations between pairs of subnetworks
% Must first run 01_mk_EEG_fMRI_corr_matrices.m

rest = 0; % 0 = use full data; 1 = use rest data only
beta_split = 1; % 0 = don't split up beta band; 1 = split up beta band
data_dir = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/correlation_data/corr_mats'];
mkdir('figs');
sublist = dir([data_dir '/sub*']);
sublist = {sublist.name};
networks = {'DNa';'DNb';'FPCNa';'FPCNb';'dATNa';'dATNb'};
%networks = {'RH_DefaultA_IPL'; 'RH_DefaultB_IPL'; 'RH_ContB_IPL'};
%networks = {'LH_DefaultA_PFCm'; 'LH_DefaultB_PFCd'; 'LH_ContB_PFCmp'};
%networks = {'LH_DefaultA_IPL'; 'LH_DefaultB_IPL'; 'LH_ContB_IPL'};
%networks = {'RH_DefaultB_Temp'; 'RH_ContB_Temp'};
if beta_split==0
    freqs = {'δ'; 'θ'; 'α'; 'β'; 'γ'};
    freqs_plot = {'δ: 1-3 Hz'; 'θ: 4-7 Hz'; 'α: 8-12: Hz'; 'β1: 13-30 Hz'; 'γ: 30-40 Hz'};
elseif beta_split==1
    freqs = {'δ'; 'θ'; 'α'; 'β1'; 'β2'; 'γ'};
    freqs_plot = {'δ: 1-3 Hz'; 'θ: 4-7 Hz'; 'α: 8-12: Hz'; 'β1: 13-20 Hz'; 'β2: 21-30 Hz'; 'γ: 31-40 Hz'};
end
network_colors = [([187 55 56]/255); ([254 147 134]/255); ([79 130 181]/255);...
    ([165 218 244]/255); ([55 119 62]/255); ([206 224 164]/255)];
do_network_pairs = 1; % set to 1 to do stats/plot network pairs
network_pairs = {[1 2];[3 4];[5 6]};
network_thr ={[-.07 .07];[-.07 .07];[-.11 .11];[-.11 .11];[-.11 .11];[-.11 .11]};
pair_names = {'DNa > DNb'; 'FPCNa > FPCNb'; 'dATNa > dATNb'};

colors = cbrewer('div','RdBu',256);
colors = abs(colors);
colors = flipud(colors);
inter_colors = cbrewer('div','PuOr',256);
inter_colors = abs(inter_colors);
inter_colors = flipud(inter_colors);

% Load in and concatenate matrices for each network
for i=1:length(networks)
    network_allsubs = []; p_mask = []; network_allsubs_lagmean = [];
    for j=1:length(sublist)
        subject_mean = [];
        if rest==1
            load([data_dir '/' sublist{j} '/' networks{i} '_mean_canonical_rest.mat']);
        elseif rest==0
            if beta_split==0
                load([data_dir '/' sublist{j} '/' networks{i} '_mean_canonical.mat']);
            elseif beta_split==1
                load([data_dir '/' sublist{j} '/' networks{i} '_mean_canonical_beta.mat']); 
                subject_mean_canonical = subject_mean_canonical_beta;
            end
        end
        subject_mean_canonical(7:end,:) = []; % chop out lags > 10 sec
        subject_mean = fisherz(subject_mean_canonical);
        network_allsubs(:,:,j) = subject_mean;
        network_allsubs_lagmean(j,:) = mean(subject_mean(1:6,:)); % mean across lags 2-8 sec
        allnetworks_allsubs{i} = network_allsubs; % store all networks
        allnetworks_allsubs_lagmean{i} = network_allsubs_lagmean;
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

% compute mean corr per network for each frequency band (mean for 0-10 sec lag)
for i=1:length(freqs)
    for j=1:length(networks)
        mean_per_freq(i,j) = mean(allnetworks_allsubs_lagmean{1,j}(:,i));
        std_per_freq(i,j) = std(allnetworks_allsubs_lagmean{1,j}(:,i));
        ste_per_freq(i,j) = std_per_freq(i,j)/(sqrt(length(allnetworks_allsubs_lagmean{1,j}(:,i))));
        lagmean_byfreq{i}(:,j) = allnetworks_allsubs_lagmean{1,j}(:,i);
    end
end

% plot mean heatmaps across all subjects for each network
    if rest==1
        figure('Position',[200,200,1500,150]);
    else
       figure('Position',[200,200,650,900]);
    end
for i=1:length(network_mean)
    if rest==1
    subplot(1,6,i)
    else
    subplot(3,2,i)
    end
    imagesc(network_mean{i}); % Display the image
    colormap(colors);
    set(gca, 'YDir', 'normal') % flip y 
    c=colorbar('Location','eastoutside');
    c.FontSize=14;
    c.Label.String=['z (\rho)' '_{EEG,fMRI}'];
    set(gca,'xtick',1:size(network_mean{i},2))
    set(gca,'ytick',1:size(network_mean{i},1))
    set(gca,'xticklabel',freqs);
    set(gca,'ytick',1:size(network_mean{i},1))
    set(gca,'yticklabel',0:2:10)
    set(gca,'TickLength',[0 0])
    set(gca,'FontSize',16,'FontName','Arial');
    xlabel({['Frequency band'];[' ']});
    if mod(i, 2) == 0
        ylabel({[' '];[' '];[' ']});
    else
        ylabel({[' '];[' '];['Lag relative to BOLD (s)']});
    end
    %c_max = max(abs(network_mean{i}(:)));
    %caxis([-c_max c_max]);
    caxis(network_thr{i})
    title(networks{i})
    box off
    hold on;
    [y_fdr,x_fdr]=find(network_p{i}==1);
    s=scatter(x_fdr,y_fdr,60,'k','*','MarkerEdgeColor',[255 215 0]/255);
    hold on;
end
print('-opengl','-r600','-dpng',['figs/heatmaps.png']);  
pause; close;

% plot network a vs b difference maps
if do_network_pairs==1
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
print('-opengl','-r600','-dpng',['figs/internetwork_heatmaps.png']); 
pause; close;
end

% plot mean corr per network (0-10 sec lag) for each frequency band
figure('Position',[200,200,1300,175]);
%figure('Position',[200,200,600,375]);
for i=1:length(freqs)
      subplot(1,length(freqs),i)
      scatter(mean_per_freq(i,:),1:length(networks),80,network_colors,'filled');
      hold on;
      errorbar(mean_per_freq(i,:),1:length(networks),ste_per_freq(i,:),'horizontal','LineWidth',1,'Color',[0 0 0],...
          'LineStyle','none');
      set(gca,'Fontsize',12,'Fontweight','normal','LineWidth',.5,'TickDir','out','box','off');
      if i==1
      yticklabels([networks]);
      else
         yticklabels([' ']); 
      end
      hold on;
      xline(0,'k--')
      xlim([-.1 .1]);
      ylim([.5 length(networks)+0.5]);
      title({[freqs_plot{i}];[' ']});
      xlabel(['z (\rho)' '_{EEG,fMRI}']);
end
print('-opengl','-r600','-dpng',['figs/freq_by_network_corr.png']); 
pause; close;
% for i=1:length(freqs)
%     subplot(1,length(freqs),i)
%     boxplot([lagmean_byfreq{1,i}],'Orientation','horizontal',...
%         'Colors',network_colors,'OutlierSize',1,'Symbol','w+',...
%         'Widths',.05,'LabelOrientation','horizontal','BoxStyle','filled','PlotStyle','compact')
%     set(gca,'Fontsize',12,'Fontweight','normal','LineWidth',.5,'TickDir','out','box','off');
%     yticklabels([networks]);
%     hold on;
%     xline(0,'k--')
% end

