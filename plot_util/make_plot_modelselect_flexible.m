function[]  = make_plot_modelselect_flexible(fcoef,coef_dim,pgon_confine,line_size)

[invx,invy] = iFD(fcoef,coef_dim);

plot(pgon_confine)
hold on 
plot(invx,invy,'r',LineWidth=line_size)

box on 
axis equal
axis([0 1 0 1])
set(gcf,'renderer','Painters')