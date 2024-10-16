function [] = make_centroid_plot_with_hist(cent,cent_fd,type)
h = scatterhist(cent_fd(:,1),cent_fd(:,2),'Kernel','on','Location','NorthEast', 'Direction','out','Marker','+',MarkerSize=5,Color='b');
hold on
if strcmp(type,'boxplot')
    hold on
    bh = boxplot(h(2),cent_fd(:,1),'orientation','horizontal','label',{''},'Widths', 0.5,'Color','b','Symbol','b+');
    set(bh,'LineWidth', 1);
%     lines = findobj(bh, 'type', 'line', 'Tag', 'Median');
%     set(lines, 'Color', 'b');

    bh = boxplot(h(3),cent_fd(:,2),'orientation','horizontal','label',{''},'Widths', 0.5,'Color','b','Symbol','b+');
    set(bh,'LineWidth', 1);
%     lines = findobj(bh, 'type', 'line', 'Tag', 'Median');
%     set(lines, 'Color', 'b');

    set(h(2:3),'XTickLabel','');
    view(h(3),[270,90]);  % Rotate the Y plot
    axis(h(1),'auto');  % Sync axes

else
    %make kernel histogram
    %h = scatterhist(cent_fd(:,1),cent_fd(:,2),'Kernel','on','Location','NorthEast', 'Direction','out','Marker','+');
    hold on
    set(h(2:3),'XTickLabel','');
    view(h(3),[270,90]);  % Rotate the Y plot
    axis(h(1),'auto');  % Sync axes
    
end
% ax = h(2);
% ax.PlotBoxAspectRatio = [1 0.2 0.2];
% ax = h(3);
% ax.PlotBoxAspectRatio = [1 0.2 0.2];
ax = h(1);
ax.Position = [0.100 0.100 0.6900 0.6900];

%scatter(cent(1),cent(2),70,'MarkerFaceColor','r','MarkerEdgeColor','r')
scatter(cent(1),cent(2),100,'MarkerEdgeColor','r','MarkerFaceColor','r')

hold off;

%plot the pdf 
