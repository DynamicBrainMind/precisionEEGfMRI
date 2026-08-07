% plot EEG-fMRI corrs for all frequency and up to 20 seconds of lag
% Must first run step1_mk_EEG_fMRI_corr_matrices.m

%% CHANGE DIRECTORIES BELOW FOR YOUR SYSTEM
data_dir = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/correlation_data/corr_mats'];
mkdir('figs/intersession');

%% settings
individual_plots=0;
sublist = {'sub-001','sub-002','sub-004','sub-006','sub-007','sub-008','sub-009','sub-010','sub-011','sub-012',...
    'sub-013','sub-014','sub-015','sub-016','sub-017','sub-018','sub-019','sub-020','sub-021','sub-023','sub-024','sub-025'};
%sublist = {'sub-001','sub-002','sub-004','sub-006','sub-007','sub-008','sub-009','sub-010','sub-011','sub-012',...
%    'sub-013','sub-014','sub-015','sub-016','sub-018','sub-019','sub-020','sub-021','sub-023','sub-024'};
networks = {'DNa';'DNb';'FPCNa';'FPCNb';'dATNa';'dATNb'};
freqs = {'δ'; 'θ'; 'α'; 'β'; 'γ'};
colors = cbrewer('div','RdBu',256);
colors = abs(colors);
colors = flipud(colors);
inter_colors = cbrewer('div','PuOr',256);
inter_colors = abs(inter_colors);
inter_colors = flipud(inter_colors);
network_colors = [([187 55 56]/255); ([254 147 134]/255); ([79 130 181]/255);...
    ([165 218 244]/255); ([55 119 62]/255); ([206 224 164]/255)]; 

%% Load in and concatenate matrices for each network
for i=1:length(networks)
    network_allsubs = []; p_mask = [];
    for j=1:length(sublist)
        subject_mean = [];
        load([data_dir '/' sublist{j} '/ses-001/' networks{i} '_mean.mat']);
        network_allsubs_ses001(:,:,j) = session_mean;
        allnetworks_allsubs_ses001{i} = network_allsubs_ses001; % store all networks
        load([data_dir '/' sublist{j} '/ses-002/' networks{i} '_mean.mat']);
        network_allsubs_ses002(:,:,j) = session_mean;
        allnetworks_allsubs_ses002{i} = network_allsubs_ses002; % store all networks
    end
end

%% compute inter-session ICC across subjects for each network
for i=1:length(networks)
    reliability_corr_matrix = NaN(size(network_allsubs_ses001,1),size(network_allsubs_ses001,2));
    reliability_ICC_matrix = NaN(size(network_allsubs_ses001,1),size(network_allsubs_ses001,2));
    % loop through each lag + frequency
    for j=1:size(network_allsubs_ses001,1) % frequencies
        for k=1:size(network_allsubs_ses001,2) % lags
            curr_sess1_vals = squeeze(allnetworks_allsubs_ses001{i}(j,k,:))';
            curr_sess2_vals = squeeze(allnetworks_allsubs_ses002{i}(j,k,:))';
            reliability_corr_matrix(j,k) = corr(curr_sess1_vals,curr_sess2_vals);
            reliability_ICC_matrix(j,k) = ICC([curr_sess1_vals,curr_sess2_vals],'C-1');
        end
    end
    allnetworks_reliability_corr_matrix{i}=reliability_corr_matrix;
    allnetworks_reliability_ICC_matrix{i}=reliability_ICC_matrix;
end

%% plot heatmaps for the two sessions for each subject and network
% make one plot for each subject
for i=1:length(sublist)
    for j=1:length(networks)
        % compute ICC
        intersess_r(i,j) = corr(fisherz(reshape(allnetworks_allsubs_ses001{j}(:,:,i),[],1)), fisherz(reshape(allnetworks_allsubs_ses002{j}(:,:,i),[],1)));
        if individual_plots==1
        % plot
        figure('Position', [100, 100, 700, 300]);
        t=tiledlayout(1,2);
        title(t,[sublist{i} ': ' networks{j} ' (r = ' num2str(intersess_r(i,j)) ')'],'FontWeight','bold','FontSize',14);
        nexttile
        imagesc(allnetworks_allsubs_ses001{j}(:,:,i)); % plot 1st session
        title(['Session 1'])
        colormap(colors);
        set(gca, 'YDir', 'normal') % flip y 
        c=colorbar('Location','southoutside');
        c.FontSize=16;
        c.Label.String=[{['\rho' '_{EEG,fMRI}'];[' ']}];
        c.Label.VerticalAlignment = 'baseline';
        set(gca,'xtick',0:5:size(allnetworks_allsubs_ses001{j}(:,:,i),2))
        set(gca,'ytick',1:size(allnetworks_allsubs_ses001{j}(:,:,i),1))
        set(gca,'ytick',1:size(allnetworks_allsubs_ses001{j}(:,:,i),1))
        set(gca,'yticklabel',0:2:20)
        set(gca,'TickLength',[0 0])
        set(gca,'FontSize',14,'FontName','Arial');
        xlabel({['Frequency'];[' ']});
        ylabel(['Lag relative to BOLD (s)']);
        %c_max = prctile(abs(allnetworks_allsubs_ses001{j}(:)),99);
        %caxis([-c_max c_max]);
        %c_limits = c.Limits;
        %c.Ticks = c_limits;
        box off
        hold on;
        nexttile
        imagesc(allnetworks_allsubs_ses002{j}(:,:,i)); % plot 2nd session
        title(['Session 2'])
        colormap(colors);
        set(gca, 'YDir', 'normal') % flip y 
        c=colorbar('Location','southoutside');
        c.FontSize=16;
        c.Label.String=[{['\rho' '_{EEG,fMRI}'];[' ']}];
        c.Label.VerticalAlignment = 'baseline';
        set(gca,'xtick',0:5:size(allnetworks_allsubs_ses001{j}(:,:,i),2))
        set(gca,'ytick',1:size(allnetworks_allsubs_ses001{j}(:,:,i),1))
        set(gca,'ytick',1:size(allnetworks_allsubs_ses001{j}(:,:,i),1))
        set(gca,'yticklabel',0:2:20)
        set(gca,'TickLength',[0 0])
        set(gca,'FontSize',14,'FontName','Arial');
        xlabel({['Frequency (Hz)'];[' ']});
        ylabel(['Lag relative to BOLD (s)']);
        %c_max = max(abs(allnetworks_allsubs_ses001{j}(:)));
        %caxis([-c_max c_max]);
        %c_limits = c.Limits;
        %c.Ticks = c_limits;
        hold on;
        print('-opengl','-r600','-dpng',['figs/intersession/' sublist{i} '_' networks{j}  '.png']); 
        close;
        end
    end
end


%% plot inter-session r values for each subject+network
data = intersess_r;
figure('Position', [100, 100, 700, 300]);
hold on;
% Prepare data for boxplot
vals = data(:);
groups = repelem(1:6, size(data,1))';
% Create boxplot (no outlier symbols)
boxplot(vals, groups, 'Colors', 'k', 'Symbol', '');
% Get box handles and fill with colors
h = findobj(gca,'Tag','Box');
colors = network_colors; % 6 distinct colors
% Boxes are returned in reverse order → flip index
for j = 1:length(h)
    patch(get(h(j),'XData'), get(h(j),'YData'), ...
          colors(end-j+1,:), ...
          'FaceAlpha', 0.75, 'EdgeColor', 'k');
end
% Overlay individual data points (jittered)
for i = 1:6
    x = i + (rand(size(data,1),1)-0.5)*0.2;
    scatter(x, data(:,i), 30, 'k', 'filled', ...
        'MarkerFaceAlpha', 0.6);
end

%% Summary metrics (mean ± SEM)
means = mean(data);
sems  = std(data) ./ sqrt(size(data,1));
% Mean markers
scatter(1:6, means, 80, 'k', 'd', 'filled');
% Error bars
errorbar(1:6, means, sems, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
% Labels and formatting
ylabel(['\itr\rm' '_{sess1,sess2}']);
title('Inter-session reliability of EEG-fMRI correlations');
xticks(1:6);
xticklabels(networks);
ylim([-1 1])
yline(0,'k-')
set(findall(gcf,'-property','FontSize'),'FontSize',18)
box on;
print('-opengl','-r600','-dpng',['figs/intersession_correlations.png']); 
pause; close;

%% plot reliability correlation matrices
colors = cbrewer('div','RdBu',256);
colors = abs(colors);
colors = flipud(colors);
tile_indices=[1 4 2 5 3 6];
figure('Position',[200,200,1200,600]);
    t = tiledlayout(2, 3);
for i=1:length(networks)
    nexttile(tile_indices(i));
    imagesc(allnetworks_reliability_ICC_matrix{i},[-.43 .86]); % Display the image
    colormap(colors);
    set(gca, 'YDir', 'normal') % flip y 
    c=colorbar('Location','eastoutside');
    c.FontSize=16;
    title(networks{i})
    c.Label.String=[{['ICC'];[' ']}];
    c.Label.VerticalAlignment = 'bottom';
    set(gca,'xtick',0:5:size(allnetworks_reliability_corr_matrix{i},2))
    set(gca,'ytick',1:size(allnetworks_reliability_corr_matrix{i},1))
    set(gca,'ytick',1:size(allnetworks_reliability_corr_matrix{i},1))
    set(gca,'yticklabel',0:2:20)
    set(gca,'TickLength',[0 0])
    set(gca,'FontSize',18,'FontName','Arial');
    xlabel({['Frequency'];[' ']});
        ylabel({[' '];[' '];['Lag relative to BOLD (s)']});

    c_limits = c.Limits;
    c.Ticks = c_limits;
    box off
    hold on;
end
print('-opengl','-r600','-dpng',['figs/ICC_maps_allnetworks.png']);  
pause; close;




