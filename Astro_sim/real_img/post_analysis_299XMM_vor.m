close all
clear
addpath(genpath('~/src/gsrg/Arp299XMM'))
addpath(genpath('~/src/gsrg/gcr_util'))
addpath(genpath('~/src/gsrg/plot_util'))
addpath(genpath('~/src/gsrg/Astro_sim/util'))
addpath(genpath('~/src/gsrg/G-SRG'))
addpath(genpath('~/src/gsrg/contour_util'))
addpath(genpath('~/src/gsrg/FD_util'))
cd '~/src/gsrg/Astro_sim/real_img/'
load('model_select_result_299XMM.mat')
load('tune_param_Arp299.mat')
load('Arp299_MOS1_evt_0.5-8.0keV_scaled.mat')
load('bootstrap_result_299XMM_vk.mat')
load('flexible_model_select_result_299XMM.mat')
%use the flexible model results
min_select_AIC = candidate_aic;
min_select_BIC = candidate_bic;
uncertain_plot_include_x = false;
uncertain_plot_include_y = false;

imagename = '~/src/gsrg/Astro_sim/real_img/sim_results/Arp299XMM_vor_vk_';
target_idx = 1;
X = double(unique(X,'rows'));
dataname = 'Arp299XMM';
keep_offset = false;

convertrate_x = 35.2;
convertrate_y = 35.2;
convertrate = convertrate_x*convertrate_y*100;
major_axis_limit = 50;
minor_axis_limit = 10;

calc_ra = @(x) calc_sec(x,7040,21640.5,25160.5,0.05,2.80783,1);
calc_dec = @(x) calc_sec(x,7044.0,24321.5,27843.5,0.05,-4.02374,2);
rightCorner = [calc_ra(-60),calc_dec(60)];
leftCorner = [calc_ra(60),calc_dec(-60)];


[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
[V, R] = voronoiDiagram(DT);


[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
[V, R] = voronoiDiagram(DT);


%% plot the histogram of voronoi area 
figure
histogram(cell_area)
saveas(gcf, strcat(imagename, 'cellarea'),'epsc')

%% scatter plot 
figure;
scatter(X(:,1),X(:,2),1,'k.')
gca = image_rescaling(gca,dataname,true,true,keep_offset);
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
gca = image_rescaling(gca,dataname,true,false,keep_offset);
saveas(gcf, strcat(imagename, 'seeds'),'epsc')
%%

%% plot segmentation result
figure;
plot_segmentation_wo_voronoi(DT, selected, cx, cy,colors ,8,true)
gca = image_rescaling(gca,dataname,true,false,keep_offset);
saveas(gcf, strcat(imagename, 'segmentation'),'epsc')

%%

%% plot contours 
line_size = 2;
figure;
glb_contour = get_contour(n, DT,selected_nonempty,adj_mat,invalid,cell_area);

target = [target_x{target_idx},target_y{target_idx}];
make_plot_contour(X,glb_contour,1,line_size,target);
gca = image_rescaling(gca,dataname,true,false,keep_offset);
saveas(gcf, strcat(imagename, 'contours'),'epsc')

%%

[centx,centy] = centroid(polyshape(target_x{target_idx},target_y{target_idx}));
[x_sample,y_sample] = sample_curve(target_x{target_idx},target_y{target_idx},300,centx, false);
fourier_coef = FD(x_sample,y_sample);

%% plot the model selection results
figure ;
make_plot_modelselect(smooth_polyset,1);
gca = image_rescaling(gca,dataname,true,true,keep_offset);
saveas(gcf, strcat(imagename, 'model_selection_AIC'),'epsc')

figure;
make_plot_modelselect(smooth_polyset,2);
gca = image_rescaling(gca,dataname,true,true,keep_offset);
saveas(gcf, strcat(imagename, 'model_selection_BIC'),'epsc')

%% plot the flexible model selection results 
figure ;
make_plot_modelselect_flexible(fourier_coef,candidate_aic,pgon_confine,line_size);
gca = image_rescaling(gca,dataname,uncertain_plot_include_x,uncertain_plot_include_y,keep_offset);
saveas(gcf, strcat(imagename, 'flexible_model_selection_AIC'),'epsc')

figure;
make_plot_modelselect_flexible(fourier_coef,candidate_bic,pgon_confine,line_size);
gca = image_rescaling(gca,dataname,uncertain_plot_include_x,uncertain_plot_include_y,keep_offset);
saveas(gcf, strcat(imagename, 'flexible_model_selection_BIC'),'epsc')

%% plot the global confidence region
r = 1024;
colors = [[0.6,0.6,0.8];[0.8,0.8,1]];
figure;
make_plot_GCR_color(X,fourier_coef,fourier_coef_all_AIC.',[1,2],min_pair_num{1}(target_idx),1000,r,colors);
gca = image_rescaling(gca,strcat(dataname,'GCR'),true,true,keep_offset);
saveas(gcf, strcat(imagename, 'GCR_AIC'),'epsc')

figure;
make_plot_GCR_color(X,fourier_coef,fourier_coef_all_BIC.',[1,2],min_pair_num{2}(target_idx),1000,r,colors);
gca = image_rescaling(gca,strcat(dataname,'GCR'),true,true,keep_offset);
saveas(gcf, strcat(imagename, 'GCR_BIC'),'epsc')

%% plot the GCR using flexible model selection
r = 1024;
colors = [[0.6,0.6,0.8];[0.8,0.8,1]];
figure;
make_plot_GCR_color_single(X,fourier_coef,fourier_coef_all_AIC.',[1,2],candidate_aic,1000,r,colors);
gca = image_rescaling(gca,strcat(dataname,'GCR'),uncertain_plot_include_x,uncertain_plot_include_y,keep_offset);
saveas(gcf, strcat(imagename, 'GCR_AIC_flexible'),'epsc')

figure;
make_plot_GCR_color_single(X,fourier_coef,fourier_coef_all_BIC.',[1,2],candidate_bic,1000,r,colors);
gca = image_rescaling(gca,strcat(dataname,'GCR'),uncertain_plot_include_x,uncertain_plot_include_y,keep_offset);
saveas(gcf, strcat(imagename, 'GCR_BIC_flexible'),'epsc')
%% plot cdf 
relfreq = true;
fontsize = 20;
legdsize = 15;
ticks_num=6;

legdvalue = {'bootstrap','bootstrap(AIC)','bootstrap(BIC)','obeserved','observerd(AIC)','observerd(BIC)'};
legdpos = 'eastoutside';
cdf_figuresize = [50 50 600 320];
%% plot flux (in bootstrap)
lwidth=2;
drop = true;
binnum = 20;
relfreq = true;
linetypes = {'-','-.','--'}; linecolors = {'k','k','k'};
num_in_poly = sum(inpolygon(X(:,1),X(:,2),target(:,1),target(:,2)));

figure

[invx,invy] = iFD(fourier_coef,min_select_AIC);
num_in_aic = sum(inpolygon(X(:,1),X(:,2),invx,invy));
[invx,invy] = iFD(fourier_coef,min_select_BIC);
num_in_bic = sum(inpolygon(X(:,1),X(:,2),invx,invy));
flux_btp{1} = num_B_poly_all_AIC;
flux_btp{2} = num_B_fd_all_AIC;
flux_btp{3} = num_B_fd_all_BIC;

make_freq_stack(flux_btp,[num_in_poly,num_in_aic,num_in_bic], ...
    legdvalue, legdpos,linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51)
set(gcf, 'Position', cdf_figuresize); %
xlabel('Flux')
ylabel('f(Flux)')
saveas(gcf, strcat(imagename, 'bootstrap_flux_AIC'),'epsc')




%% plot flux (in observe)
drop = true;
binnum = 20;
relfreq = true;
linetypes = {'-','-.','--'}; linecolors = {'k','k','k'};
num_in_poly = sum(inpolygon(X(:,1),X(:,2),target(:,1),target(:,2)));

figure

[invx,invy] = iFD(fourier_coef,min_select_AIC);
num_in_aic = sum(inpolygon(X(:,1),X(:,2),invx,invy));
[invx,invy] = iFD(fourier_coef,min_select_BIC);
num_in_bic = sum(inpolygon(X(:,1),X(:,2),invx,invy));
flux_btp{1} = num_O_poly_all_AIC;
flux_btp{2} = num_O_fd_all_AIC;
flux_btp{3} = num_O_fd_all_BIC;

xlabel('Flux')
ylabel('f(Flux)')
make_freq_stack(flux_btp,[num_in_poly,num_in_aic,num_in_bic], ...
    legdvalue, legdpos,linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51)
set(gcf, 'Position', cdf_figuresize); %
saveas(gcf, strcat(imagename, 'bootstrap_flux_2'),'epsc')



%% plot adjusted flux (in bootstrap)
drop = true;
binnum = 20;
relfreq = true;
linetypes = {'-','-.','--'}; linecolors = {'k','k','k'};
area_poly = polyarea(target(:,1),target(:,2));

num_in_poly = sum(inpolygon(X(:,1),X(:,2),target(:,1),target(:,2)));
num_in_poly = get_adjust_flux(num_in_poly,area_poly,n);

figure

[invx,invy] = iFD(fourier_coef,min_select_AIC);
num_in_aic= sum(inpolygon(X(:,1),X(:,2),invx,invy));
area_aic = polyarea(invx,invy);
num_in_aic = get_adjust_flux(num_in_aic,area_aic,n);
[invx,invy] = iFD(fourier_coef,min_select_BIC);
num_in_bic = sum(inpolygon(X(:,1),X(:,2),invx,invy));
area_bic = polyarea(invx,invy);
num_in_bic = get_adjust_flux(num_in_bic,area_bic,n);

flux_btp{1} = get_adjust_flux(num_B_poly_all_AIC,area_poly_all_AIC,n);
flux_btp{2} = get_adjust_flux(num_B_fd_all_AIC,area_fd_all_AIC,n);
flux_btp{3} = get_adjust_flux(num_B_fd_all_BIC,area_fd_all_BIC,n);

make_freq_stack(flux_btp,[num_in_poly,num_in_aic,num_in_bic], ...
    legdvalue, legdpos,linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51)
xlabel('Adjust Flux')
ylabel('f(Adjust Flux)')
set(gcf, 'Position', cdf_figuresize); %
saveas(gcf, strcat(imagename, 'bootstrap_adjustflux'),'epsc')


figure
cdf_figuresize_simple = [50 50 450 320];
flux_btp = {};
flux_btp{1} = get_adjust_flux(num_B_poly_all_AIC,area_poly_all_AIC,n);

make_freq_stack(flux_btp,[num_in_poly], ...
    {}, legdpos,linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51)
xlabel('Adjust Flux')
ylabel('f(Adjust Flux)')
set(gcf, 'Position', cdf_figuresize_simple); %
saveas(gcf, strcat(imagename, 'bootstrap_adjustflux_simple'),'epsc')



%% plot adjust flux (in observe)
drop = true;
binnum = 20;
relfreq = true;
linetypes = {'-','-.','--'}; linecolors = {'k','k','k'};
area_poly = polyarea(target(:,1),target(:,2));

num_in_poly = sum(inpolygon(X(:,1),X(:,2),target(:,1),target(:,2)));
num_in_poly = get_adjust_flux(num_in_poly,area_poly,n);

figure

[invx,invy] = iFD(fourier_coef   ,min_select_AIC);
num_in_aic= sum(inpolygon(X(:,1),X(:,2),invx,invy));
area_aic = polyarea(invx,invy);
num_in_aic = get_adjust_flux(num_in_aic,area_aic,n);
[invx,invy] = iFD(fourier_coef,min_select_BIC);
num_in_bic = sum(inpolygon(X(:,1),X(:,2),invx,invy));
area_bic = polyarea(invx,invy);
num_in_bic = get_adjust_flux(num_in_bic,area_bic,n);

flux_btp{1} = get_adjust_flux(num_O_poly_all_AIC,area_poly_all_AIC,n);
flux_btp{2} = get_adjust_flux(num_O_fd_all_AIC,area_fd_all_AIC,n);
flux_btp{3} = get_adjust_flux(num_O_fd_all_BIC,area_fd_all_BIC,n);

make_freq_stack(flux_btp,[num_in_poly,num_in_aic,num_in_bic], ...
    legdvalue, legdpos,linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51)
xlabel('Adjust Flux')
ylabel('f(Adjust Flux)')
set(gcf, 'Position', cdf_figuresize); %
saveas(gcf, strcat(imagename, 'bootstrap_adjustflux_2'),'epsc')





%% plot area 
drop = true;
binnum = 20;
%bwidth = 0.0008;
%xlimit = linspace(1000/convertrate,12000/convertrate,ticks_num);
area_poly = polyarea(target(:,1),target(:,2))*convertrate;
[invx,invy] = iFD(fourier_coef,min_select_AIC);
area_aic = polyarea(invx,invy)*convertrate;
[invx,invy] = iFD(fourier_coef,min_select_BIC);
area_bic = polyarea(invx,invy)*convertrate;
area_btp{1} = area_poly_all_AIC*convertrate;
area_btp{2} = area_fd_all_AIC*convertrate;
area_btp{3} = area_fd_all_BIC*convertrate;

figure
make_freq_stack(area_btp,[area_poly,area_aic,area_bic], ...
    legdvalue, legdpos,linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51)
% xticks(xlimit)
% xticklabels(round(xlimit*24*24*100))
xlabel('Area')
ylabel('f(Area)')
set(gcf, 'Position', cdf_figuresize); %
saveas(gcf, strcat(imagename, 'bootstrap_area'),'epsc')


figure
area_btp = {};
area_btp{1} = area_poly_all_AIC*convertrate;

make_freq_stack(area_btp,[area_poly], ...
    {}, legdpos,linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,51)
%xticks(xlimit)
%xticklabels(round(xlimit*24*24*100))
xlabel('Area')
ylabel('f(Area)')
set(gcf, 'Position', cdf_figuresize_simple); %
saveas(gcf, strcat(imagename, 'bootstrap_area_simple'),'epsc')



%% plot the centroid bootstrap
figure
make_centroid_plot_with_hist([centx,centy],[centx_fd_all_AIC;centy_fd_all_AIC]','boxplot')
gca = image_rescaling(gca,strcat(dataname,'centroid'));
saveas(gcf, strcat(imagename, 'centroid_scatter_box'),'epsc')

figure
make_centroid_plot_with_hist([centx,centy],[centx_fd_all_AIC;centy_fd_all_AIC]','histogram')
gca = image_rescaling(gca,strcat(dataname,'centroid'));
saveas(gcf, strcat(imagename, 'centroid_scatter_hist'),'epsc')

%%

%% plot the eccentricity 
%rescale to real data scale

fourier_coef_all_AIC_rescale = rescale_coef(fourier_coef_all_AIC.',0,0,convertrate_x,convertrate_y);
long_all_axis = abs(abs(fourier_coef_all_AIC_rescale(2,:))+abs(fourier_coef_all_AIC_rescale(end,:))); 
short_all_axis = abs(abs(fourier_coef_all_AIC_rescale(2,:))-abs(fourier_coef_all_AIC_rescale(end,:))); 
fourier_coef_rescale = rescale_coef(fourier_coef,0,0,convertrate_x,convertrate_y);
long_axis = abs(abs(fourier_coef_rescale(2))+abs(fourier_coef_rescale(end))); 
short_axis = abs(abs(fourier_coef_rescale(2))-abs(fourier_coef_rescale(end))); 
figure
make_eccentricity_plot_with_hist([long_axis,short_axis],[long_all_axis;short_all_axis]',minor_axis_limit,major_axis_limit, dataname,'boxplot')
saveas(gcf, strcat(imagename, 'eccentricity_scatter_box'),'epsc')
figure
make_eccentricity_plot_with_hist([long_axis,short_axis],[long_all_axis;short_all_axis]',minor_axis_limit,major_axis_limit, dataname,'histogram')
saveas(gcf, strcat(imagename, 'eccentricity_scatter_hist'),'epsc')



%% plot 10 illustration of bootstrap results
figure 
[centx,centy] = centroid(polyshape(target_x{target_idx},target_y{target_idx}));
[x_sample,y_sample] = sample_curve(target_x{target_idx},target_y{target_idx},300,centx, false);
plot(x_sample,y_sample,'b',LineWidth=2)
hold on
[invx,invy] = iFD(fourier_coef,min_select_AIC);
plot(invx,invy,'r',LineWidth=2)
gca = image_rescaling(gca,dataname,true,true,keep_offset);
saveas(gcf, strcat(imagename, 'illustration_fd_AIC'),'epsc')

figure 
hold on
for ii = 1:10
    [invx,invy] = iFD(fourier_coef_all_AIC(ii,:).',300);
    plot(invx,invy,'Color',GRAY,LineWidth=2)
end
plot(x_sample,y_sample,'b',LineWidth=2)
box on 
axis equal
gca = image_rescaling(gca,dataname,uncertain_plot_include_x,uncertain_plot_include_y,keep_offset);
saveas(gcf, strcat(imagename, 'illustration_bootstrap_10'),'epsc')

figure
for ii = 1:10
    subplot(2,5,ii)
    [invx,invy] = iFD(fourier_coef_all_AIC(ii,:).',300);
    hold on 
    plot(invx,invy,'Color',GRAY,LineWidth=2)
    plot(x_sample,y_sample,'b',LineWidth=2)

    gca = image_rescaling(gca,dataname,true,true,keep_offset);

end

%% plot 10 illustration of bootstrap results for BIC
figure 
[centx,centy] = centroid(polyshape(target_x{target_idx},target_y{target_idx}));
[x_sample,y_sample] = sample_curve(target_x{target_idx},target_y{target_idx},300,centx, false);
plot(x_sample,y_sample,'b',LineWidth=2)
hold on
[invx,invy] = iFD(fourier_coef,min_select_BIC);
plot(invx,invy,'r',LineWidth=2)
gca = image_rescaling(gca,dataname,true,true,keep_offset);
saveas(gcf, strcat(imagename, 'illustration_fd_BIC'),'epsc')

figure 
hold on
for ii = 1:10
    [invx,invy] = iFD(fourier_coef_all_BIC(ii,:).',300);
    plot(invx,invy,'Color',GRAY,LineWidth=2)
end

plot(x_sample,y_sample,'b',LineWidth=2)

gca = image_rescaling(gca,dataname,true,true,keep_offset);
saveas(gcf, strcat(imagename, 'illustration_bootstrap_BIC'),'epsc')

plot_kde(X,1,1,1)
gca = image_rescaling(gca,dataname,true,true,keep_offset);
cb=colorbar('eastoutside');
cb.Position = cb.Position + [0.5*1e-1,0,0,0];

saveas(gcf, strcat(imagename, 'kde_contour'),'epsc')




figure
[invx,invy] = iFD(fourier_coef,300);

plot(invx,invy,'Color',GRAY,LineWidth=2)
hold on
plot(x_sample,y_sample,'b',LineWidth=2)
figure
[invx,invy] = iFD(reshape(fourier_coef, [], 1),300);

plot(invx,invy,'Color',GRAY,LineWidth=2)
hold on
plot(x_sample,y_sample,'b',LineWidth=2)

