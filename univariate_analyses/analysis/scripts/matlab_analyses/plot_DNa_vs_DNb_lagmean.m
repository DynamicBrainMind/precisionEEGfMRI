% Plot summaries of EEG-fMRI correlations for Yeo7 DN vs DNa and DNb

data_dir = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/correlation_data/corr_mats'];
mkdir('figs');
sublist = dir([data_dir '/sub*']);
sublist = {sublist.name};
networks = {'DNa';'DNb'};
freqs = {'Delta'; 'Theta'; 'Alpha'; 'Beta1'; 'Beta2'; 'Gamma'};
network_colors = [[[187 55 56]/255]; [[254 147 134]/255]];

% Load in and concatenate matrices for each network
for i=1:length(networks)
    network_allsubs = []; p_mask = [];
    for j=1:length(sublist)
        subject_mean = [];
        load([data_dir '/' sublist{j} '/' networks{i} '_mean_canonical_beta.mat']);
        subject_mean_canonical_beta(7:end,:) = []; % chop out lags > 10 sec
        subject_mean = fisherz(subject_mean_canonical_beta);
        network_allsubs(:,:,j) = subject_mean;
        allnetworks_allsubs{i} = network_allsubs; % store all networks
    end
    % get correlations for each network at each frequency (mean across lags)
    for j=1:length(freqs)
            freq_network_corrs = [];
            freq_network_corrs = squeeze(allnetworks_allsubs{i}(:,j,:)); % all lags
            freq_network_corrs = mean(freq_network_corrs)'; % mean across lags
            all_freq_network_corrs{j}(:,i) = freq_network_corrs;     
    end
end

% two-way RM-ANOVA
Y = []; S = []; F1=[]; F2=[]; 
FACTNAMES = {'frequency';'network'};
for i=1:length(all_freq_network_corrs)
   Y = [Y; all_freq_network_corrs{i}(:)];
   S = [S; repmat(1:size(all_freq_network_corrs{i},1),1,2)'];
   F1 = [F1; ones(length(all_freq_network_corrs{i}(:)),1)*i]; % freq
   F2 = [F2; ones(size(all_freq_network_corrs{i},1),1); ones(size(all_freq_network_corrs{i},1),1)*2]; % network
end
rm_anova_stats = rm_anova2(Y,S,F1,F2,FACTNAMES);

% Stats for each frequency band
for i=1:length(freqs)
    % Friedman test for interaction
    [p_freq(i),friedman_table{i},friedman_stats{i}] = friedman(all_freq_network_corrs{i},1,'off');
    % pairwise Wilcoxon sign rank tests
    for j=1:length(networks)
        for k=1:length(networks)
        p_pairwise{i}(j,k)=signrank(all_freq_network_corrs{i}(:,j),all_freq_network_corrs{i}(:,k));
        end
        p_pairwise{i}(eye(size(p_pairwise{i})) == 1) = NaN; 
        upper_triangle_with_zeros = triu(p_pairwise{i});
        lower_triangle_mask = (upper_triangle_with_zeros == 0);
        p_pairwise{i}(lower_triangle_mask) = NaN;
    end
    % apply FDR correction
    p_fdr{i} = p_pairwise{i}(:);
    p_fdr{i}(isnan(p_fdr{i}))=[]; % remove NaNs so that fdr ignores them
    p_fdr{i} = fdr(p_fdr{i});
end 

% organize data for RM-ANOVA



% Get mean + STE for plotting
for i=1:length(freqs)
    mean_corrs{i} = mean(all_freq_network_corrs{i});
    n = size(all_freq_network_corrs{i}, 1); % Number of observations (rows)
    std_M = std(all_freq_network_corrs{i}, 0, 1); % Standard deviation for each column
    ste_corrs{i} = std_M / sqrt(n); % Standard error for each column
end

%% Plot
figure('Position',[200,200,600,600]);
x_start = 1;
for i=1:length(freqs)
    scatter([x_start:x_start+size(all_freq_network_corrs{i},2)-1],mean_corrs{i},90,network_colors,'filled')
    hold on;
    errorbar([x_start:x_start+size(all_freq_network_corrs{i},2)-1],mean_corrs{i},ste_corrs{i},'LineWidth',1,'Color',[0 0 0]);
    set(gca,'Fontsize',22,'Fontweight','normal','LineWidth',1,'TickDir','out','box','off');
    yline(0,'k--');
    %xlabel('Frequency band','FontSize',26);
    ylabel(['z (\rho)' '_{EEG,fMRI}'],'FontSize',26)
hold on;
    x_starts(i) = x_start+1;
    x_start = x_start + size(all_freq_network_corrs{i},2)+1;
    hold on
end
hold on
set(gca,'Xtick',x_starts-.5,'XTickLabel',{'δ [1-3 Hz]','θ [4-7 Hz]', 'α [8-12 Hz]','β1 [13-20 Hz]','β2 [20-30 Hz]','γ [31-40 Hz]'});
xtickangle(45)
% significant pairwise effects 
hold on;
xcoords = [[1,2];[1,3];[2,3]]; % needs editing if not using 3 networks
for i=1:length(p_fdr)
    coords = x_starts(i)+xcoords-2;
    for j=1:length(p_fdr{i})
        if p_fdr{i}(j)<0.05
        sigstar(coords(j,:), p_fdr{i}(j),0,[0 0 0]); 
    end
    end
end
ylim([-0.03 0.12])
print('-opengl','-r600','-dpng',['figs/DNa_vs_DNb_MSHBM_eeg_fmri.png']);  


