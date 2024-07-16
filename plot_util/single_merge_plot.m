
function ax = single_merge_plot(ax, metrics, gridsizes,qt,color, offset,LineStyle,Marker,linename)

K = length(gridsizes);
metrics = reshape(metrics,K,[]);
medval =  median(metrics,2)';

%plot((1:K)+offset,medval,'LineStyle',LineStyle,'Marker',Marker,'Color',color)
errorbar((1:K)+offset,medval,prctile(metrics,qt(1),2)'-medval,prctile(metrics,qt(2),2)'-medval,'LineStyle',LineStyle,'Marker',Marker,'Color',color,'DisplayName',sprintf(linename))
% for i =1:K
%     boxplot(metrics(i,:), 'positions', i+offset(i),'Color',color);
% end