function [] = make_hist_stack(data, meanval,edges,leg_title,legendpos)
map = brewermap(3,'Set1'); 
%figure
hold on
%for i = 1:length(data)
%histogram(data{i},edges,'facecolor',map(i,:),'facealpha',.5,'edgecolor','none','Normalization','probability')

%end
histogram(data{1},edges,'facecolor',map(1,:),'facealpha',.5,'edgecolor','none','Normalization','probability')
histogram(data{2},edges,'facecolor',map(3,:),'facealpha',.5,'edgecolor','none','Normalization','probability')
yl = ylim;
xline(meanval(1), 'Color', 'b', 'LineWidth', 2);
xline(meanval(2), 'Color', 'm', 'LineWidth', 2);
box on
axis tight
legend(leg_title,'Location',legendpos)
%legalpha('H1','H2','location','northwest')