function plot_eeg_fmri_individuals(network_allsubs,network_name)

colors = cbrewer('div','RdBu',256);
colors = abs(colors);
colors = flipud(colors);
freqs = {'δ'; 'θ'; 'α'; 'β1'; 'β2'; 'γ'};

figure('Position', [100, 100, 1200, 1200]);
tiledlayout(5,5)
for i=1:size(network_allsubs,3)
    nexttile
    imagesc(network_allsubs(:,:,i)); % Display the image
    colormap(colors);
    set(gca, 'YDir', 'normal') % flip y 
    c=colorbar('Location','eastoutside');
    c.FontSize=10;
    c.Label.String=['z (\rho)' '_{EEG,fMRI}'];
    set(gca,'xtick',1:size(network_allsubs(:,:,i),2))
    set(gca,'ytick',1:size(network_allsubs(:,:,i),1))
    set(gca,'xticklabel',freqs);
    set(gca,'ytick',1:size(network_allsubs(:,:,i),1))
    set(gca,'yticklabel',0:2:10)
    set(gca,'TickLength',[0 0])
    set(gca,'FontSize',10,'FontName','Arial');
    xlabel({['Frequency band'];[' ']});
    %if mod(i, 2) == 0
    %    ylabel({[' '];[' '];[' ']});
    %else
        ylabel({[' '];[' '];['Lag relative to BOLD (s)']});
    %end
    c_max = max(abs(reshape(network_allsubs(:,:,i),[],1)));
    caxis([-c_max c_max]);
    %caxis(network_thr{i})
    title(network_name)
    box off
    %hold on;
    %[y_fdr,x_fdr]=find(network_p{i}==1);
    %s=scatter(x_fdr,y_fdr,60,'k','*','MarkerEdgeColor',[255 215 0]/255);
    hold on;
end
pause; close;