%% Generate time/drug fraction occupancy profiles for BRAFi/MEKi
%% Miles Miller, Thomas Ng

close all

namer = 'DAB_TRAM'
tail  = '';

brafname = 'Dab';
mekname  = 'Tram';

BRAF_TSUM = DAB_tfrac;
MEK_TSUM  = TRAM_tfrac;

%Normalize to maximum frac. occ
BRAF_norm = BRAF_TSUM./max(BRAF_TSUM(:));
MEK_norm  = MEK_TSUM./max(MEK_TSUM(:));

map = colormap(jet);
C = 1:size(BRAF_norm,1);
C = (C-1);
C = C./max(C);

%% Plot profile for tumor core
figure, arrowPlot_CC(mean(BRAF_norm(:,1:10),2), mean(MEK_norm(:,1:10),2), C', map)
xlabel(['[' brafname '] Frac. occ']), ylabel(['[' mekname '] Frac. occ']), 
xlim([-0.05 1.05]), ylim([-0.05 1.05])
title('Core');

xticks([0 0.5 1]);
yticks([0 0.5 1]);
colormap(map)
d = colorbar;
d.Ticks = [0 0.5 1];
caxis([min(C), max(C)]);
set(gca,'linewidth',1.5)
saveas(gcf, fullfile('',[namer '_core_' tail '.pdf']));

%% Plot profile for tumor edge
figure,  arrowPlot_CC(mean(BRAF_norm(:,40:50),2), mean(MEK_norm(:,40:50),2), C',map)
xlabel(['[' brafname '] Frac. occ']), ylabel(['[' mekname '] Frac. occ']),
xlim([-0.05 1.05]), ylim([-0.05 1.05])
title('Edge')

xticks([0 0.5 1]);
yticks([0 0.5 1]);
colormap(map)
d = colorbar;
d.Ticks = [0 0.5 1];
caxis([min(C), max(C)]);
set(gca,'linewidth',1.5)
saveas(gcf, fullfile('',[namer '_edge_' tail '.pdf']));

return