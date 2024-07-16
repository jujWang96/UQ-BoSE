
close all
clear
addpath(genpath('/Users/joycewang/Desktop/gsrg/Arp299'))
addpath(genpath('~/Desktop/gsrg/gcr_util'))
addpath(genpath('~/Desktop/gsrg/plot_util'))
addpath(genpath('~/Desktop/gsrg/Astro_sim/util'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))
addpath(genpath('~/Desktop/gsrg/contour_util'))
addpath(genpath('~/Desktop/gsrg/FD_util'))
cd '~'/Desktop/gsrg/Astro_sim/real_img/
load('model_select_result_299XMM.mat')
load('tune_param_Arp299.mat')
load('Arp299_MOS1_evt_0.5-8.0keV_scaled.mat')
load('bootstrap_result_299XMM.mat')

imagename = 'Arp299XMM_';
target_idx = 1;
X = double(unique(X,'rows'));


calc_ra = @(x) calc_sec(x,7040,21640.5,25160.5,0.05,2.80783,1);
calc_dec = @(x) calc_sec(x,7044.0,24321.5,27843.5,0.05,-4.02374,2);
rightCorner = [calc_ra(-150),calc_dec(150)];
leftCorner = [calc_ra(150),calc_dec(-150)];
fov = polyshape([leftCorner(1),rightCorner(1),rightCorner(1),leftCorner(1)], ...
    [leftCorner(2),leftCorner(2),rightCorner(2),rightCorner(2)]);
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
[V, R] = voronoiDiagram(DT);

%% scatter plot 
figure;
scatter(X(:,1),X(:,2),1,'k.')
gca = image_rescaling(gca,'Arp299XMM');
saveas(gcf, strcat(imagename, 'scatterplot'),'epsc')
disp(strcat('the total flux in fov is ',num2str(sum(inpolygon(X(:,1),X(:,2), ...
    [leftCorner(1),rightCorner(1),rightCorner(1),leftCorner(1)],[leftCorner(2),leftCorner(2),rightCorner(2),rightCorner(2)])))))

%%


% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
%[seeds, num] = get_seeds_sim_kde(cx, cy,invalid, adj_mat, threshold);
[seeds, num] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat,cell_area,P, threshold);

%% plot seeds
figure;
plot_seeds(DT, cx, cy, seeds, [], [], lines(num), num, 0)
gca = image_rescaling(gca,'Arp299XMM');
saveas(gcf, strcat(imagename, 'seeds'),'epsc')
%%
%% plot segmentation result
figure;
plot_segmentation_wo_voronoi(DT, selected, cx, cy,colors ,8,true)

gca = image_rescaling(gca,'Arp299XMM');
saveas(gcf, strcat(imagename, 'segmentation'),'epsc')

%% plot contours 
figure;
glb_contour = get_contour(n, DT,selected_nonempty,adj_mat,invalid,cell_area);

target = [target_x{target_idx},target_y{target_idx}];
make_plot_contour(X,glb_contour,1,2,target);
gca = image_rescaling(gca,'Arp299XMM');
saveas(gcf, strcat(imagename, 'contours'),'epsc')
%%

[centx,centy] = centroid(polyshape(target_x{target_idx},target_y{target_idx}));
[x_sample,y_sample] = sample_curve(target_x{target_idx},target_y{target_idx},300,centx, false);
fourier_coef = FD(x_sample,y_sample)';

%% plot the model selection results
figure ;
make_plot_modelselect(smooth_polyset,1);
gca = image_rescaling(gca,'Arp299XMM');
saveas(gcf, strcat(imagename, 'model_selection_AIC'),'epsc')

figure;
make_plot_modelselect(smooth_polyset,2);
gca = image_rescaling(gca,'Arp299XMM');
saveas(gcf, strcat(imagename, 'model_selection_BIC'),'epsc')


%% plot the global confidence region
r = 1024;
colors = [[0.6,0.6,0.8];[0.8,0.8,1]];
figure;
make_plot_GCR_color(X,fourier_coef',fourier_coef_all_AIC',[1,2],min_pair_num{1}(1),1000,r,colors);
gca = image_rescaling(gca,'Arp299XMMGCR');
saveas(gcf, strcat(imagename, 'GCR_AIC'),'epsc')

figure;
make_plot_GCR_color(X,fourier_coef',fourier_coef_all_BIC',[1,2],min_pair_num{2}(1),1000,r,colors);
gca = image_rescaling(gca,'Arp299XMMGCR');
saveas(gcf, strcat(imagename, 'GCR_BIC'),'epsc')





%% plot flux (in bootstrap)
drop = true;
relfreq = true;
fontsize = 15;
legdsize = 20;
binnum = 30;
lwidth = 1.5;
linetypes = {'-','-.'}; linecolors = {'k','k'};
num_in_poly = sum(inpolygon(X(:,1),X(:,2),target(:,1),target(:,2)));

figure
min_select_AIC
[invx,invy] = iFD(fourier_coef',min_select_AIC);
num_in_fd = sum(inpolygon(X(:,1),X(:,2),invx,invy));
flux_aic{1} = num_B_poly_all_AIC;
flux_aic{2} = num_B_fd_all_AIC;

make_freq_stack(flux_aic,[num_in_poly,num_in_fd], ...
    {'bootstrap flux','bootstrap flux (AIC)','obeserved flux','observerd AIC flux'}, ...
    'northwest',linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51)
set(gcf, 'Position', [50 50 600 400]); %
saveas(gcf, strcat(imagename, 'bootstrap_flux_AIC'),'epsc')

figure
make_cdf_stack(flux_aic,[num_in_poly,num_in_fd], ...
    {'bootstrap value','bootstrap AIC value','obeserved value','observerd AIC value'}, ...
    'northwest',linetypes,linecolors,20,lwidth,fontsize,drop,legdsize)
%xlim([2500,3000])
xlabel('Flux')
ylabel('F(Flux)')
set(gcf, 'Position', [50 50 600 400]); 
saveas(gcf, strcat(imagename, 'bootstrap_flux_cdf_AIC'),'epsc')
fprintf("The 5 percentile of bootstrap flux %f\n",prctile(flux_aic{1},5))
fprintf("The 95 percentile of bootstrap flux %f\n",prctile(flux_aic{1},95))
fprintf("The 5 percentile of bootstrap aic flux %f\n",prctile(flux_aic{2},5))
fprintf("The 95 percentile of bootstrap aic flux %f\n",prctile(flux_aic{2},95))


figure
min_select_BIC
[invx,invy] = iFD(fourier_coef',min_select_BIC);
num_in_fd = sum(inpolygon(X(:,1),X(:,2),invx,invy));
flux_bic{1} = num_B_poly_all_BIC;
flux_bic{2} = num_B_fd_all_BIC;

make_freq_stack( flux_bic,[num_in_poly,num_in_fd], ...
    {'bootstrap flux','bootstrap flux (BIC)','obeserved flux','observerd BIC flux'}, ...
    'northwest',linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51);

set(gcf, 'Position', [50 50 600 400]); %
saveas(gcf, strcat(imagename, 'bootstrap_flux_BIC'),'epsc')
figure
make_cdf_stack(flux_bic,[num_in_poly,num_in_fd], ...
    {'bootstrap value','bootstrap BIC value','obeserved value','observerd BIC value'}, ...
    'northwest',linetypes,linecolors,20,lwidth,fontsize,drop,legdsize)
%xlim([2500,3000])
xlabel('Flux')
ylabel('F(Flux)')
set(gcf, 'Position', [50 50 600 400]); 
saveas(gcf, strcat(imagename, 'bootstrap_flux_cdf_BIC'),'epsc')



%% plot adjusted flux (in bootstrap)
drop = true;
relfreq = true;
fontsize = 15;
legdsize = 20;
binnum = 30;
lwidth = 1.5;
linetypes = {'-','-.'}; linecolors = {'k','k'};
area_poly = polyarea(target(:,1),target(:,2));
num_in_poly = sum(inpolygon(X(:,1),X(:,2),target(:,1),target(:,2)));

num_in_poly = get_adjust_flux(num_in_poly,area_poly,n);
figure
min_select_AIC

num_in_fd = sum(inpolygon(X(:,1),X(:,2),invx,invy));
[invx,invy] = iFD(fourier_coef',min_select_AIC);
area_fd = polyarea(invx,invy);
num_in_fd = get_adjust_flux(num_in_fd,area_fd,n);

flux_aic{1} = get_adjust_flux(num_B_poly_all_AIC,area_poly_all_AIC,n);
flux_aic{2} = get_adjust_flux(num_B_fd_all_AIC,area_fd_all_AIC,n);

make_freq_stack(flux_aic,[num_in_poly,num_in_fd], ...
    {'bootstrap flux','bootstrap flux (AIC)','obeserved flux','observerd AIC flux'}, ...
    'northwest',linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51)
set(gcf, 'Position', [50 50 600 400]); %
saveas(gcf, strcat(imagename, 'bootstrap_adjustflux_AIC'),'epsc')

figure
make_cdf_stack(flux_aic,[num_in_poly,num_in_fd], ...
    {'bootstrap value','bootstrap AIC value','obeserved value','observerd AIC value'}, ...
    'northwest',linetypes,linecolors,20,lwidth,fontsize,drop,legdsize)
%xlim([2500,3000])
xlabel('Flux')
ylabel('F(Flux)')
set(gcf, 'Position', [50 50 600 400]); 
saveas(gcf, strcat(imagename, 'bootstrap_adjustflux_cdf_AIC'),'epsc')
fprintf("The 5 percentile of bootstrap flux %f\n",prctile(flux_aic{1},5))
fprintf("The 95 percentile of bootstrap flux %f\n",prctile(flux_aic{1},95))
fprintf("The 5 percentile of bootstrap aic flux %f\n",prctile(flux_aic{2},5))
fprintf("The 95 percentile of bootstrap aic flux %f\n",prctile(flux_aic{2},95))


figure
min_select_BIC
[invx,invy] = iFD(fourier_coef',min_select_BIC);
num_in_fd = sum(inpolygon(X(:,1),X(:,2),invx,invy));
area_fd = polyarea(invx,invy);
num_in_fd = get_adjust_flux(num_in_fd,area_fd,n);
flux_aic{1} = get_adjust_flux(num_B_poly_all_BIC,area_poly_all_BIC,n);
flux_aic{2} = get_adjust_flux(num_B_fd_all_BIC,area_fd_all_BIC,n);

make_freq_stack( flux_bic,[num_in_poly,num_in_fd], ...
    {'bootstrap flux','bootstrap flux (BIC)','obeserved flux','observerd BIC flux'}, ...
    'northwest',linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51);

set(gcf, 'Position', [50 50 600 400]); %
saveas(gcf, strcat(imagename, 'bootstrap_adjustflux_BIC'),'epsc')
figure
make_cdf_stack(flux_bic,[num_in_poly,num_in_fd], ...
    {'bootstrap value','bootstrap BIC value','obeserved value','observerd BIC value'}, ...
    'northwest',linetypes,linecolors,20,lwidth,fontsize,drop,legdsize)
%xlim([2500,3000])
xlabel('Flux')
ylabel('F(Flux)')
set(gcf, 'Position', [50 50 600 400]); 
saveas(gcf, strcat(imagename, 'bootstrap_adjustflux_cdf_BIC'),'epsc')



%%

%% plot flux (in observe)
drop = true;
relfreq = true;
fontsize = 15;
legdsize = 20;
binnum = 30;
lwidth = 1.5;
linetypes = {'-','-.'}; linecolors = {'k','k'};
num_in_poly = sum(inpolygon(X(:,1),X(:,2),target(:,1),target(:,2)));

figure
min_select_AIC
[invx,invy] = iFD(fourier_coef',min_select_AIC);
num_in_fd = sum(inpolygon(X(:,1),X(:,2),invx,invy));
flux_aic{1} = num_O_poly_all_AIC;
flux_aic{2} = num_O_fd_all_AIC;

make_freq_stack(flux_aic,[num_in_poly,num_in_fd], ...
    {'bootstrap flux','bootstrap flux (AIC)','obeserved flux','observerd AIC flux'}, ...
    'northwest',linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51)
set(gcf, 'Position', [50 50 600 400]); %
saveas(gcf, strcat(imagename, 'bootstrap_flux_AIC2'),'epsc')

figure
make_cdf_stack(flux_aic,[num_in_poly,num_in_fd], ...
    {'bootstrap value','bootstrap AIC value','obeserved value','observerd AIC value'}, ...
    'northwest',linetypes,linecolors,20,lwidth,fontsize,drop,legdsize)
%xlim([2500,3000])
xlabel('Flux')
ylabel('F(Flux)')
set(gcf, 'Position', [50 50 600 400]); 
saveas(gcf, strcat(imagename, 'bootstrap_flux_cdf_AIC2'),'epsc')
fprintf("The 5 percentile of bootstrap flux %f\n",prctile(flux_aic{1},5))
fprintf("The 95 percentile of bootstrap flux %f\n",prctile(flux_aic{1},95))
fprintf("The 5 percentile of bootstrap aic flux %f\n",prctile(flux_aic{2},5))
fprintf("The 95 percentile of bootstrap aic flux %f\n",prctile(flux_aic{2},95))


figure
min_select_BIC
[invx,invy] = iFD(fourier_coef',min_select_BIC);
num_in_fd = sum(inpolygon(X(:,1),X(:,2),invx,invy));
flux_bic{1} = num_O_poly_all_BIC;
flux_bic{2} = num_O_fd_all_BIC;

make_freq_stack( flux_bic,[num_in_poly,num_in_fd], ...
    {'bootstrap flux','bootstrap flux (BIC)','obeserved flux','observerd BIC flux'}, ...
    'northwest',linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51);

set(gcf, 'Position', [50 50 600 400]); %
saveas(gcf, strcat(imagename, 'bootstrap_flux_BIC2'),'epsc')
figure
make_cdf_stack(flux_bic,[num_in_poly,num_in_fd], ...
    {'bootstrap value','bootstrap BIC value','obeserved value','observerd BIC value'}, ...
    'northwest',linetypes,linecolors,20,lwidth,fontsize,drop,legdsize)
set(gcf, 'Position', [50 50 600 400]); 
%xlim([2500,3000])
xlabel('Flux')
ylabel('F(Flux)')
saveas(gcf, strcat(imagename, 'bootstrap_flux_cdf_BIC2'),'epsc')



%% plot area 
bwidth = 0.001;
drop = true;
relfreq = true;
binnum = 30;
convertrate = 35.2*35.2*100;
xlimit = linspace(1000/convertrate,4700/convertrate,21);

area_poly = polyarea(target(:,1),target(:,2));
[invx,invy] = iFD(fourier_coef',min_select_AIC);
area_fd = polyarea(invx,invy);
area_aic{1} = area_poly_all_AIC;
area_aic{2} = area_fd_all_AIC;
figure
make_freq_stack(area_aic,[area_poly,area_fd], ...
    {'bootstrap area','bootstrap area (AIC)','obeserved area','observerd AIC area'}, ...
    'northeast',linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51)
xticks(xlimit)
xticklabels(round(xlimit*convertrate))
set(gcf, 'Position', [50 50 600 400]); %
saveas(gcf, strcat(imagename, 'bootstrap_area_AIC'),'epsc')

figure
make_cdf_stack(area_aic,[area_poly,area_fd], ...
    {'bootstrap value','bootstrap AIC value','obeserved value','observerd AIC value'}, ...
    'northwest',linetypes,linecolors,bwidth,lwidth,fontsize,drop,legdsize)
xlimit = [3900,5000];
xticks(linspace(xlimit(1),xlimit(2),8)/convertrate)
xticklabels(round(linspace(xlimit(1),xlimit(2),8)))

xlim(xlimit/convertrate)
xlabel('Area')
ylabel('F(Area)')
set(gcf, 'Position', [50 50 600 400]); %
saveas(gcf, strcat(imagename, 'bootstrap_area_cdf_AIC'),'epsc')


[invx,invy] = iFD(fourier_coef',min_select_BIC);
area_fd = polyarea(invx,invy);
area_bic{1} = area_poly_all_BIC;
area_bic{2} = area_fd_all_BIC;
figure
make_freq_stack(area_bic,[area_poly,area_fd], ...
    {'bootstrap area','bootstrap area (BIC)','obeserved area','observerd BIC area'}, ...
    'northwest',linetypes,linecolors,binnum,lwidth, fontsize,drop,legdsize,relfreq,51)
xticks(xlimit)
xticklabels(round(xlimit*convertrate))
set(gcf, 'Position', [50 50 600 400]); %
saveas(gcf, strcat(imagename, 'bootstrap_area_BIC'),'epsc')

figure
make_cdf_stack(area_bic,[area_poly,area_fd], ...
    {'bootstrap value','bootstrap BIC value','obeserved value','observerd BIC value'}, ...
    'southeast',linetypes,linecolors,bwidth,lwidth, fontsize,drop,legdsize)
xlimit = [3900,5000];
xticks(linspace(xlimit(1),xlimit(2),8)/convertrate)
xticklabels(round(linspace(xlimit(1),xlimit(2),8)))

xlim(xlimit/convertrate)
xlabel('Area')
ylabel('F(Area)')
set(gcf, 'Position', [50 50 600 400]); %
saveas(gcf, strcat(imagename, 'bootstrap_area_cdf_BIC'),'epsc')


%% plot the centroid bootstrap
figure
make_centroid_plot([centx,centy],[centx_fd_all_AIC;abs(centy_fd_all_AIC)]',0.5)
gca = image_rescaling(gca,'Arp299XMM');
saveas(gcf, strcat(imagename, 'centroid_scatter_AIC'),'epsc')

%%

%% plot the eccentricity 
%rescale to real data scale
figure
fourier_coef_all_AIC_rescale = rescale_coef(fourier_coef_all_AIC',0,0,35.2,35.2);
long_all_axis = abs(abs(fourier_coef_all_AIC_rescale(2,:))+abs(fourier_coef_all_AIC_rescale(end,:))); 
short_all_axis = abs(abs(fourier_coef_all_AIC_rescale(2,:))-abs(fourier_coef_all_AIC_rescale(end,:))); 
fourier_coef_rescale = rescale_coef(fourier_coef,0,0,35.2,35.2);
long_axis = abs(abs(fourier_coef_rescale(2))+abs(fourier_coef_rescale(end))); 
short_axis = abs(abs(fourier_coef_rescale(2))-abs(fourier_coef_rescale(end))); 
make_eccentricity_plot([long_axis,short_axis],[long_all_axis;short_all_axis]',30,45, 'Arp299XMM',0.5)
saveas(gcf, strcat(imagename, 'eccentricity_scatter_AIC'),'epsc')

%%


%% plot the centroid bootstrap
figure
make_centroid_plot([centx,centy],[centx_fd_all_BIC;abs(centy_fd_all_BIC)]',0.5)
gca = image_rescaling(gca,'Arp299XMM');
saveas(gcf, strcat(imagename, 'centroid_scatter_BIC'),'epsc')

%%

%% plot the eccentricity 
%rescale to real data scale
figure
fourier_coef_all_BIC_rescale = rescale_coef(fourier_coef_all_BIC',0,0, 35.2, 35.2);
long_all_axis = abs(abs(fourier_coef_all_BIC_rescale(2,:))+abs(fourier_coef_all_BIC_rescale(end,:))); 
short_all_axis = abs(abs(fourier_coef_all_BIC_rescale(2,:))-abs(fourier_coef_all_BIC_rescale(end,:))); 
fourier_coef_rescale = rescale_coef(fourier_coef,0,0, 35.2, 35.2);
long_axis = abs(abs(fourier_coef_rescale(2))+abs(fourier_coef_rescale(end))); 
short_axis = abs(abs(fourier_coef_rescale(2))-abs(fourier_coef_rescale(end))); 
make_eccentricity_plot([long_axis,short_axis],[long_all_axis;short_all_axis]',30,45, 'Arp299XMM',0.5)
saveas(gcf, strcat(imagename, 'eccentricity_scatter_BIC'),'epsc')

%% plot 10 illustration of bootstrap results
figure 
[centx,centy] = centroid(polyshape(target_x{target_idx},target_y{target_idx}));
[x_sample,y_sample] = sample_curve(target_x{target_idx},target_y{target_idx},300,centx, false);
plot(x_sample,y_sample,'b',LineWidth=2)
hold on
[invx,invy] = iFD(fourier_coef',min_select_AIC);
plot(invx,invy,'r',LineWidth=2)
gca = image_rescaling(gca,'Arp299XMM');
saveas(gcf, strcat(imagename, 'illustration_fd_AIC'),'epsc')

figure 
hold on
for ii = 1:10
    [invx,invy] = iFD(fourier_coef_all_AIC(ii,:).',300);
    plot(invx,invy,'Color',GRAY,LineWidth=2)
end
plot(x_sample,y_sample,'b',LineWidth=2)

gca = image_rescaling(gca,'Arp299XMM');
saveas(gcf, strcat(imagename, 'illustration_bootstrap_AIC'),'epsc')


%% plot 10 illustration of bootstrap results for BIC
figure 
[centx,centy] = centroid(polyshape(target_x{target_idx},target_y{target_idx}));
[x_sample,y_sample] = sample_curve(target_x{target_idx},target_y{target_idx},300,centx, false);
plot(x_sample,y_sample,'b',LineWidth=2)
hold on
[invx,invy] = iFD(fourier_coef',min_select_BIC);
plot(invx,invy,'r',LineWidth=2)
gca = image_rescaling(gca,'Arp299XMM');
saveas(gcf, strcat(imagename, 'illustration_fd_BIC'),'epsc')

figure 
hold on
for ii = 1:10
    [invx,invy] = iFD(fourier_coef_all_BIC(ii,:).',300);
    plot(invx,invy,'Color',GRAY,LineWidth=2)
end
plot(x_sample,y_sample,'b',LineWidth=2)

gca = image_rescaling(gca,'Arp299XMM');
saveas(gcf, strcat(imagename, 'illustration_bootstrap_BIC'),'epsc')









