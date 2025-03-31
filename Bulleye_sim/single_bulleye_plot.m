function ax = single_bulleye_plot(ax,x,linestyle, linecolor,displacename)
%pd = fitdist(x,'Kernel');
%xgrid = linspace(0,0.3,100)';
%pdfEst = pdf(pd,xgrid);
%line(xgrid,pdfEst)

[f,xi] = ksdensity(x,'Support','positive','BoundaryCorrection','reflection','Bandwidth',0.01);
plot(xi,f,LineStyle=linestyle,Color=linecolor,DisplayName=displacename);
xlim([0,0.3])
box on

