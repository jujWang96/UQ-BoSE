function [] = make_freq_stack( data, meanval,leg_title,legendpos,linestyles,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,axisnum,xlimit)
hold on
maxl = 0;
if drop
    %drop the top 1% and bottom 1% values of distribution
    for ii = 1:length(data)
        upb = prctile(data{ii},99);
        lowb = prctile(data{ii},1);
        data{ii} = data{ii}(data{ii}>lowb & data{ii}<upb);
    end

end
for ii = 1:length(data)
    maxl = max(maxl,length(data{ii}));

    h = histfit(data{ii}',binnum,'Kernel');
    %delete(h(1))
    h(1).FaceColor = "#80B3FF";
    h(1).FaceAlpha = 0.5;
    set(h(2),'color',linecolors{ii},'LineStyle',linestyles{ii},'LineWidth', 1)
end
for ii = 1:length(meanval)
    xline(meanval(ii), 'Color', linecolors{ii}, 'LineWidth', lwidth,'Color','r');
end
box on
axis tight

if relfreq
    %yticklabels(linspace(0,maxl,axisnum))
    
    N1=sum(h(1).YData);
    h(1).YData=h(1).YData/N1;
    h(2).YData=h(2).YData/N1;
    yticks(linspace(0,1,axisnum))
%     yticklabels(linspace(0,1,axisnum))
end
if length(leg_title)>1
    legend(leg_title,'Location',legendpos,'FontSize',legdsize)
end
set(gca,'FontSize',fontsize)
if nargin == 14
    xlim(xlimit)
end
