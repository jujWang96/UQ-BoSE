function [] = make_centroid_plot(cent,cent_fd,alpha)

scatter(cent_fd(:,1),cent_fd(:,2),'MarkerFaceColor','b','MarkerFaceAlpha',alpha,'MarkerEdgeAlpha',alpha)
hold on
scatter(cent(1),cent(2),70,'MarkerFaceColor','r','MarkerEdgeColor','r')

%plot the pdf 
