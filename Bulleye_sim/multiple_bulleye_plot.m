function ax = multiple_bulleye_plot(ax,x,linestyles,linecolors,displacenames,addtitle, titlename,addxlabel,labelname)
hold on
for i = 1:length(x)
    ax = single_bulleye_plot(ax,x{i},linestyles{i},linecolors{i},displacenames{i});
end
if addtitle
    title(titlename)
    %set(gca,'xaxisLocation','top')
end
if addxlabel
    xlabel(labelname)
end