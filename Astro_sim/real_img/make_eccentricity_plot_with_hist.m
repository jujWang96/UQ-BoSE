function [] = make_eccentricity_plot_with_hist(ecce,ecce_fd,llim,ulim, inputname,type)

h = scatterhist(ecce_fd(:,1),ecce_fd(:,2),'Kernel','on','Location','NorthEast', 'Direction','out','Marker','+',MarkerSize=5,Color='b');
hold on
if strcmp(type,'boxplot')
    hold on
    bh = boxplot(h(2),ecce_fd(:,1),'orientation','horizontal','label',{''},'Widths', 0.5,'Color','b','Symbol','b+');
    set(bh,'LineWidth', 1);

    bh = boxplot(h(3),ecce_fd(:,2),'orientation','horizontal','label',{''},'Widths', 0.5,'Color','b','Symbol','b+');
    set(bh,'LineWidth', 1);

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
ax = h(1);
ax.Position = [0.100 0.100 0.6900 0.6900];
scatter(ecce(1),ecce(2),100,'MarkerFaceColor','r','MarkerEdgeColor', 'r')
axis([llim ulim llim ulim ])
box on
axis square
if strcmp('NGC2300',inputname)
    text(0.35*(ulim-llim)+llim, 0.9*(ulim-llim)+llim, 'NGC2300 (Chandra)','FontSize',18,'FontWeight','bold')
elseif strcmp('NGC2300XMM',inputname)
    text(0.41*(ulim-llim)+llim, 0.9*(ulim-llim)+llim, 'NGC2300 (XMM)','FontSize',18,'FontWeight','bold')
else
    text(0.5*(ulim-llim)+llim, 0.9*(ulim-llim)+llim, 'ARP299 (XMM)','FontSize',18,'FontWeight','bold')

end
if strcmp('NGC2300',inputname)

    Ylm=ylim;                          % get x, y axis limits 
    Xlm=xlim;                          % so can position relative instead of absolute
    Xlb=mean(Xlm);                    % set horizontally at midpoint
    Ylb=0.99*Ylm(1);      
    xlabel('Major semiaxis [arcsec]', ...
        'Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center')
    Xlb=0.9*Xlm(1);                    % set horizontally at midpoint
    Ylb = mean(Ylm); 
    ylabel('Minor semiaxis [arcsec]','Position',[Xlb Ylb], ...
        'VerticalAlignment','baseline','HorizontalAlignment','center')
else

    Ylm=ylim;                          % get x, y axis limits 
    Xlm=xlim;                          % so can position relative instead of absolute
    Xlb=mean(Xlm);                    % set horizontally at midpoint
    Ylb=0.99*Ylm(1);      
    xlabel('Major semiaxis [arcsec]', ...
        'Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center')
    Xlb=0.9*Xlm(1);                    % set horizontally at midpoint
    Ylb = mean(Ylm); 
    ylabel('Minor semiaxis [arcsec]','Position',[Xlb Ylb], ...
        'VerticalAlignment','baseline','HorizontalAlignment','center')

end
set(gca,'XTickLabelRotation', 0,'linewidth',2)
set(gca,'Fontsize',18,'FontWeight','bold')
set(gcf, 'Position', [50 50 400 400]); %
