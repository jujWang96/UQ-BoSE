function [] = make_dist_stack( data, meanval,leg_title,legendpos,linestyles,linecolors,bwidth,lwidth,fontsize,drop,legdsize)

hold on
if drop
    %drop the top 1% and bottom 1% values of distribution
    upb = prctile(data{1},99);
    lowb = prctile(data{1},1);
    data{1} = data{1}(data{1}>lowb & data{1}<upb);
    upb = prctile(data{2},99);
    lowb = prctile(data{2},1);
    data{2} = data{2}(data{2}>lowb & data{2}<upb);

end

[f,xi] = ksdensity(data{1},'Support','positive', ...
    'BoundaryCorrection','reflection','Bandwidth',bwidth);
plot(xi,f,LineStyle=linestyles{1},Color=linecolors{1},LineWidth=lwidth);

[f,xi] = ksdensity(data{2},'Support','positive', ...
    'BoundaryCorrection','reflection','Bandwidth',bwidth);
plot(xi,f,LineStyle=linestyles{2},Color=linecolors{2},LineWidth=lwidth);
xline(meanval(1), 'Color', 'b', 'LineWidth', lwidth);
xline(meanval(2), 'Color', 'm', 'LineWidth', lwidth);
box on
axis tight
legend(leg_title,'Location',legendpos,'FontSize',legdsize)
set(gca,'FontSize',fontsize)
