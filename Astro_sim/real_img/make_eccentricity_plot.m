function [] = make_eccentricity_plot(ecce,ecce_fd,llim,ulim, inputname,alpha)

scatter(ecce_fd(:,1),ecce_fd(:,2),'MarkerFaceColor','b','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
hold on
scatter(ecce(1),ecce(2),70,'MarkerFaceColor','r','MarkerEdgeColor', 'r')
axis([llim ulim llim ulim ])
box on
axis square
if strcmp('NGC2300',inputname)
    text(0.72*(ulim-llim)+llim, 0.92*(ulim-llim)+llim, inputname,'FontSize',18,'FontWeight','bold')
elseif strcmp('NGC2300XMM',inputname)
    text(0.57*(ulim-llim)+llim, 0.92*(ulim-llim)+llim, inputname,'FontSize',18,'FontWeight','bold')
else
    text(0.57*(ulim-llim)+llim, 0.92*(ulim-llim)+llim, inputname,'FontSize',18,'FontWeight','bold')

end
set(gca,'XTickLabelRotation', 0,'linewidth',2)
set(gca,'Fontsize',18,'FontWeight','bold')
set(gcf, 'Position', [50 50 400 400]); %
