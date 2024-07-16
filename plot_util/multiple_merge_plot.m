
function ax = multiple_merge_plot(ax, metrics,metrics_beam,metrics_random, gridsizes,qt,colors, offsets,LineStyle,Marker,linenames)
hold on 
ax = single_merge_plot(ax, metrics, gridsizes,qt,colors(1),offsets(1,:),LineStyle,Marker,linenames{1});
ax = single_merge_plot(ax, metrics_beam, gridsizes,qt,colors(2),offsets(2,:),LineStyle,Marker,linenames{2});
ax = single_merge_plot(ax, metrics_random, gridsizes,qt,colors(3),offsets(3,:),LineStyle,Marker,linenames{3});
K = length(gridsizes);
xticks(1:K)
xticklabels(gridsizes)
xlim([1+min(offsets(:))-0.1,K+max(offsets(:))+0.1])
box on