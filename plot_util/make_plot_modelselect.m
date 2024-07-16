function[]  = make_plot_modelselect(polyset,idx)
pgon = polyset{idx};
hold on 
for i =1:length(pgon)
    plot(pgon{i})
end
box on 
axis equal
axis([0 1 0 1])
set(gcf,'renderer','Painters')
