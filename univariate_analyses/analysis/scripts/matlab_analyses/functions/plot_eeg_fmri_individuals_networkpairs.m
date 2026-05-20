function plot_eeg_fmri_individuals(allnetworks_allsubs,networks,sublist)

mkdir('figs/individual')
colors = cbrewer('div','RdBu',256);
colors = abs(colors);
colors = flipud(colors);
freqs = {'δ'; 'θ'; 'α'; 'β1'; 'β2'; 'γ'};
network_thr ={[-.07 .07];[-.07 .07];[-.11 .11];[-.11 .11];[-.11 .11];[-.11 .11]};

for i=1:size(allnetworks_allsubs{1},3) % for each subject
    for j=1:size(allnetworks_allsubs,2)% for each network pair
    if mod(j, 2) ~= 0 % plot from odd networks only to set up pairs
        figure('Position', [100, 100, 700, 300]);
        t=tiledlayout(1,2)
        title(t,sublist{i},'FontWeight','bold','FontSize',14);
        %t.Title.FontWeight = 'bold';
    nexttile
        imagesc(allnetworks_allsubs{j}(:,:,i)); % plot 1st network pair
        title(networks{j})
        colormap(colors);
        set(gca, 'YDir', 'normal') % flip y 
        c=colorbar('Location','southoutside');
        c.FontSize=10;
        c.Label.String=['z (\rho)' '_{EEG,fMRI}'];
        set(gca,'xtick',1:size(allnetworks_allsubs{j}(:,:,i),2))
        set(gca,'ytick',1:size(allnetworks_allsubs{j}(:,:,i),1))
        set(gca,'xticklabel',freqs);
        set(gca,'ytick',1:size(allnetworks_allsubs{j}(:,:,i),1))
        set(gca,'yticklabel',0:2:10)
        set(gca,'TickLength',[0 0])
        set(gca,'FontSize',14,'FontName','Arial');
        xlabel({['Frequency band'];[' ']});
        ylabel({[' '];[' '];['Lag relative to BOLD (s)']});
        caxis(network_thr{j})
        box off
        hold on;
    nexttile
        imagesc(allnetworks_allsubs{j+1}(:,:,i)); % plot 2nd network pair
        title(networks{j+1})
        colormap(colors);
        set(gca, 'YDir', 'normal') % flip y 
        c=colorbar('Location','southoutside');
        c.FontSize=10;
        c.Label.String=[['z (\rho)' '_{EEG,fMRI}']];
        set(gca,'xtick',1:size(allnetworks_allsubs{j}(:,:,i),2))
        set(gca,'ytick',1:size(allnetworks_allsubs{j}(:,:,i),1))
        set(gca,'xticklabel',freqs);
        set(gca,'ytick',1:size(allnetworks_allsubs{j}(:,:,i),1))
        set(gca,'yticklabel',0:2:10)
        set(gca,'TickLength',[0 0])
        set(gca,'FontSize',14,'FontName','Arial');
        xlabel({['Frequency band'];[' ']});
        ylabel({[' '];[' '];['Lag relative to BOLD (s)']});
        caxis(network_thr{j})
        set(gcf, 'PaperPositionMode', 'auto');
        print('-opengl','-r600','-dpng',['figs/individual/' sublist{i} '_' networks{j}  '.png']); 
    close;
    end
    
    end

end