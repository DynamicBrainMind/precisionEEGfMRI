% Plot summaries of ICC for each network as a function of lag
% x axis = lag; y axis = ICC

%% load ICC for each network (for every channel)
data_dir = ['/Users/ak4379/Documents/data/R21_EEG-fMRI/derivatives/correlation_data'];
data = importdata([data_dir '/iccs_full.csv']);
ICC_all = [NaN; data.data];
textdata = data.textdata;
networks = {'DNa';'DNb';'FPCNa';'FPCNb';'dATNa';'dATNb'};
lag_secs = {'0';'2';'4';'6';'8';'10';'12';'14';'16';'18';'20'};
stackbar_colors = [[1 1 1]; [.5 .5 .5]; [0 0 0]];
icc_excellent = 0.75;
icc_good = 0.6;
icc_fair = 0.4;
icc_poor = 0;

%% Get all ICCs for each network at each lag
lags = unique(data.textdata(2:end,3));
lags = sort(str2double(lags));
for i=1:length(networks)
    network_inds = strmatch(networks{i},textdata(:,4),'exact');
    ICC_networks{i} = [];
    for j=1:length(lags)
        curr_lag = lags(j);
        lag_inds = strmatch(num2str(curr_lag),textdata(:,3),'exact');
        network_lag_inds = intersect(network_inds,lag_inds);
        ICC_networks{i} = [ICC_networks{i} ICC_all(network_lag_inds)];
    end
end

%% find % ICC excellent, good, fair, poor per lag (per network)
for i=1:length(ICC_networks)
    for j=1:length(lags)
        excellent_perc_network{j}(i) = length(find(ICC_networks{i}(:,j)>icc_excellent))/size(ICC_networks{i}(:,j),1);
        good_perc_network{j}(i) = length(find(ICC_networks{i}(:,j)<icc_excellent & ICC_networks{i}(:,j)>icc_good))/size(ICC_networks{i}(:,j),1);
        fair_perc_network{j}(i) = length(find(ICC_networks{i}(:,j)<icc_good & ICC_networks{i}(:,j)>icc_fair))/size(ICC_networks{i}(:,j),1);
        poor_perc_network{j}(i) = length(find(ICC_networks{i}(:,j)<icc_fair))/size(ICC_networks{i}(:,j),1);
    end
end

%% find overall % ICC excellent, good, fair, poor per lag
% combine ICCs across networks for each lag
for i=1:length(lags)
    ICC_combined{i} = [];
    for j=1:length(networks)
        ICC_combined{i} = [ICC_combined{i}; ICC_networks{j}(:,i)];
    end
    excellent_perc_overall(i) = length(find(ICC_combined{i}>icc_excellent))/size(ICC_combined{i},1);
    good_perc_overall(i) = length(find(ICC_combined{i}<icc_excellent & ICC_combined{i}>icc_good))/size(ICC_combined{i},1);
    fair_perc_overall(i) = length(find(ICC_combined{i}<icc_good & ICC_combined{i}>icc_fair))/size(ICC_combined{i},1);
    poor_perc_overall(i) = length(find(ICC_combined{i}<icc_fair))/size(ICC_combined{i},1);
end

% organize data for stacked bar plot
stacked_icc = [good_perc_overall; fair_perc_overall; excellent_perc_overall]';

%% Plot ICC by lag for each network
figure('Position',[200,200,600,250]);
h = bar(stacked_icc,'stacked')
set(gca,'Fontsize',16,'Fontweight','normal','LineWidth',1,'TickDir','out','box','off');
ylabel('Proportion of ICC values');
xlabel('Lag (sec) relative to BOLD');
xticklabels(lag_secs);
ylim([0 .6]);
for i = 1:numel(h)
    h(i).FaceColor = stackbar_colors(i,:);
end
legend('Fair ICC (0.40-0.59)','Good ICC (0.60-0.74)','Excellent ICC (0.75-1.00)')
print('-opengl','-r600','-dpng',['figs/icc_lag.png']);  

