%% Generate time/drug profiles for BRAFi/MEKi
%% Miles Miller, Thomas Ng

close all

namer = 'VEM_COBI'              % Filename
tail  = '';                     % suffix

brafname = 'Vem';
mekname  = 'Cobi';

BRAF_TSUM = VEM_tsum;
MEK_TSUM  = COBI_tsum;

%Normalize to maximum [drug]
BRAF_norm = BRAF_TSUM./max(BRAF_TSUM(:));
MEK_norm  = MEK_TSUM./max(MEK_TSUM(:));

map = colormap(jet);
C = 1:size(BRAF_TSUM,1);
C = (C-1);
C = C./max(C);

%% Plot profile for tumor core
figure, arrowPlot_CC(mean(BRAF_norm(:,1:10),2), mean(MEK_norm(:,1:10),2), C', map)
xlabel([brafname ' Frac. max']), ylabel([mekname ' Frac. max']), 
xlim([-0.05 1.05]), ylim([-0.05 1.05])

xticks([0 0.5 1]);
yticks([0 0.5 1]);
colormap(map)
d = colorbar;
d.Ticks = [0 0.5 1];
d.TickLabels = {'0', '' , '15 days'};

caxis([min(C), max(C)]);
set(gca,'linewidth',1.5)
title('Core')
saveas(gcf, fullfile('',[namer '_core_' tail '.pdf']));

%% Plot profile for tumor edge
figure,  arrowPlot_CC(mean(BRAF_norm(:,40:50),2), mean(MEK_norm(:,40:50),2), C',map)
xlabel([brafname ' Frac. max']), ylabel([mekname ' Frac. max']), 
xlim([-0.05 1.05]), ylim([-0.05 1.05])

xticks([0 0.5 1]);
yticks([0 0.5 1]);
colormap(map)
d = colorbar;
d.Ticks = [0 0.5 1];
d.TickLabels = {'0', '' , '15 days'};
caxis([min(C), max(C)]);
set(gca,'linewidth',1.5)
title('Edge')
saveas(gcf, fullfile('',[namer '_edge_' tail '.pdf']));


return
