% run null models
% must also separately run compute_eegfmri_correlations.m (for real data) to generate plots

n_perm = 1000;
n_bands = 6;
networks = {'DNa'; 'DNb'};
mkdir('figs');

signed_network_diff_all = zeros(n_bands,n_perm);
abs_network_diff_all = zeros(n_bands,n_perm);

parpool(6);
parfor i = 1:n_perm
   perm = num2str(i);
    r_sub_avg_DNa = compute_eegfmri_correlations_null('DNa',perm);
    r_sub_avg_DNb = compute_eegfmri_correlations_null('DNb',perm);
    signed_network_diff = r_sub_avg_DNa-r_sub_avg_DNb;
    abs_network_diff = abs(signed_network_diff);
    signed_network_diff_all(:,i) = mean(signed_network_diff,2);
    abs_network_diff_all(:,i) = mean(abs_network_diff,2);
    fprintf('done permutation %d\n',i);
end

% save permutation results
mkdir('results');
save(['results/null_distribution_abs.mat'],'abs_network_diff_all');
save(['results/null_distribution_signed.mat'],'signed_network_diff_all');

% compare against real data
load(['results/real_distribution.mat']);
real_diffs_abs = mean(abs_network_diff,2);
load(['results/real_distribution_signed.mat']);
real_diffs_signed = mean(network_diff,2);

% compute p for each frequency (signed)
for i = 1:length(real_diffs_signed)
    all_diffs = []; sorted_diffs = [];
    all_diffs=[real_diffs_signed(i); signed_network_diff_all(i,:)'];
    sorted_diffs=sort(all_diffs,'descend');
    position_true(i)=find(sorted_diffs==real_diffs_signed(i));
    p_vals_signed(i)=position_true(i)/n_perm;
end

% plot null vs real for each frequency (signed)
freq_names = {'δ (1-3 Hz)'; 'θ (4-7 Hz)'; 'α (8-12 Hz)'; 'β1 (13-20 Hz)'; 'β2 (21-30 Hz)'; 'γ (31-40 Hz)'};
figure('Position',[200,200,1000,1000]);
tiledlayout(3,2)
for i = 1:length(real_diffs_signed)
    nexttile
    %histogram(signed_network_diff_all(i,:), 'Normalization', 'pdf', 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
    histogram(signed_network_diff_all(i,:), 'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none');
    hold on;
    xline(real_diffs_signed(i), 'Color', [.7, .13, .13], 'LineWidth', 2, 'Label', ['Observed \Delta\rho'],'FontSize',14);
    xlabel({['Surrogate inter-network difference in']; ['EEG-BOLD correlation (\Delta\rho)']});
    ylabel('Number of Permutations');
    title([freq_names{i}]);
    set(gca,'Fontsize',16,'Fontweight','normal','LineWidth',.5,'TickDir','out','box','off');
    %xlim([-.01 real_diffs(i)+.01]);
    xlim([-.045 .075]);
    hold on;
end
print('-opengl','-r300','-dpng',['figs/null_models.png']);


