
close all
clear
addpath(genpath('~/src/gsrg/NGC2300'))
addpath(genpath('~/src/gsrg/gcr_util'))
addpath(genpath('~/src/gsrg/plot_util'))
addpath(genpath('~/src/gsrg/Astro_sim/util'))
addpath(genpath('~/src/gsrg/G-SRG'))
addpath(genpath('~/src/gsrg/contour_util'))
addpath(genpath('~/src/gsrg/FD_util'))
cd '~/src/gsrg/Astro_sim/real_img/'
load('tune_param_2300chandra.mat')
load('~/src/gsrg/NGC2300/ngc2300_box_058kev_evt2.mat')

load('model_select_result_2300chandra.mat')

load('bootstrap_result_2300chandra_th15_vk.mat')
load('flexible_model_select_result_2300chandra.mat')
%use the flexible model results
min_select_AIC = candidate_aic;
min_select_BIC = candidate_bic;
imagename = '~/src/gsrg/Astro_sim/real_img/sim_results/NGC2300chandra_vor_vk_';
dataname = 'NGC2300';
keep_offset = false;
uncertain_plot_include_x = false;
uncertain_plot_include_y = false;
X = double(unique(X,'rows'));

target_idx = 4;

calc_ra = @(x) calc_sec(x,706.0,3062.0,3415.0,0.492,0,1);
calc_dec = @(x) calc_sec(x,682.0,3489.0,3830.0,0.492,0,2);

rightCorner = [calc_ra(-60),calc_dec(60)];
leftCorner = [calc_ra(60),calc_dec(-60)];
fov = polyshape([leftCorner(1),rightCorner(1),rightCorner(1),leftCorner(1)], ...
    [leftCorner(2),leftCorner(2),rightCorner(2),rightCorner(2)]);

convertrate_x = 34.7;
convertrate_y = 33.5;
convertrate = convertrate_x*convertrate_y*100;
major_axis_limit = 50;
minor_axis_limit = 10;


% %% get the SNR
% area_poly = polyarea(target(:,1),target(:,2));
% num_in_poly = sum(inpolygon(X(:,1),X(:,2),target(:,1),target(:,2)));
% mu_background = (n-num_in_poly)/(1-area_poly);
% 
% num_in_poly = get_adjust_flux(num_in_poly,area_poly,n);
% mu_in = num_in_poly/area_poly;
% SNR = mu_in/sqrt(mu_in+area_poly*mu_background)
% %% get contrast
% target2 = [target_x{3},target_y{3}];
% area_poly_c = polyarea(target2(:,1),target2(:,2));
% num_in_poly_c = sum(inpolygon(X(:,1),X(:,2),target2(:,1),target2(:,2)));
% mu_in_c = num_in_poly_c/area_poly_c-mu_background;
% constast = num_in_poly_c/num_in_poly_c*mu_in_c/num_in_poly
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
figure
colors = lines(length(selected));

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
pdf_figuresize = [50 50 600 320];
%% plot flux (in bootstrap)
lwidth=1.5;
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
set(gcf, 'Position', pdf_figuresize); %
xlabel('Flux')
ylabel('f(Flux)')
saveas(gcf, strcat(imagename, 'bootstrap_flux_AIC'),'epsc')




%% plot flux (in observe)
drop = true;
lwidth = 1.5;
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
set(gcf, 'Position', pdf_figuresize); %
saveas(gcf, strcat(imagename, 'bootstrap_flux_2'),'epsc')



%% plot adjusted flux (in bootstrap)
drop = true;
lwidth = 1.5;
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
set(gcf, 'Position', pdf_figuresize); %
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
lwidth = 1.5;
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
set(gcf, 'Position', pdf_figuresize); %
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
set(gcf, 'Position', pdf_figuresize); %
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



