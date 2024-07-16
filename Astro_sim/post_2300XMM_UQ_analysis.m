clear 
close all
parpath = '~/src/gsrg/';
resultpath = 'Astro_sim/real_data/results/';
dataname = "NGC2300XMM";

load(strcat(parpath,resultpath,dataname,"_vor_UQ_analysis.mat"))
plotpath = 'Astro_sim/real_data/plots/';

addpath(genpath(strcat(parpath,'Astro_sim')))
addpath(genpath(strcat(parpath,'G-SRG')))
addpath(genpath(strcat(parpath,'plot_util')))
addpath(genpath(strcat(parpath,'contour_util')))
addpath(genpath(strcat(parpath,'FD_util')))
addpath(genpath(strcat(parpath,'boot_util')))
addpath(genpath(strcat(parpath,'gcr_util')))
addpath(genpath(strcat(parpath,'NGC2300')))

imagename = strcat(parpath, plotpath, dataname, "_vor_UQ_");


%% plot parameters
line_size = 2;

convertrate_x = 24;
convertrate_y = 24;
convertrate = convertrate_x*convertrate_y*100;
major_axis_limit = 60;
minor_axis_limit = 10;
keep_offset = false;
uncertain_plot_include_x = true;
uncertain_plot_include_y = false;

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

% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
%[seeds, num] = get_seeds_sim_kde(cx, cy,invalid, adj_mat, threshold);
[seeds, num] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat,cell_area,best_P, best_threshold);

disp(strcat("Num of seeds is ",num2str(num)))
disp(strcat("Value of FD (BIC) is ",num2str(candidate_bic)))
disp(strcat("Number of segments is ",num2str(length(raw_contour.regionId))))

%% plot seeds
figure;
plot_seeds(DT, cx, cy, seeds, [], [], lines(num), num, 0)
gca = image_rescaling(gca,dataname,true,false,keep_offset);
saveas(gcf, strcat(imagename, 'seeds'),'epsc')
%% plot segmentation result
figure;
plot_segmentation_wo_voronoi(DT, selected, cx, cy,lines(length(selected)) ,8,true);
gca = image_rescaling(gca,dataname,true,false,keep_offset);
saveas(gcf, strcat(imagename, 'segmentation'),'epsc')


%% plot contours 
figure;

target = [contourx,contoury];
make_plot_contour(X,raw_contour,1,line_size,target);
gca = image_rescaling(gca,dataname,true,false,keep_offset);
saveas(gcf, strcat(imagename, 'contours'),'epsc')

%% plot the flexible model selection results 

figure ;
make_plot_modelselect_flexible(fourier_coef,candidate_aic,pgon_confine,line_size);
gca = image_rescaling(gca,dataname,uncertain_plot_include_x,uncertain_plot_include_y,keep_offset);
saveas(gcf, strcat(imagename, 'flexible_model_selection_AIC'),'epsc')

figure;
make_plot_modelselect_flexible(fourier_coef,candidate_bic,pgon_confine,line_size);
gca = image_rescaling(gca,dataname,uncertain_plot_include_x,true,keep_offset);
saveas(gcf, strcat(imagename, 'flexible_model_selection_BIC'),'epsc')


%% plot the GCR using flexible model selection
colors = [[0.6,0.6,0.8];[0.8,0.8,1]];
figure;
make_plot_GCR_color_polygon(X,fourier_coef,fourier_coef_b_subset,[1,2],candidate_aic,1000,colors);
gca = image_rescaling(gca,dataname,uncertain_plot_include_x,uncertain_plot_include_y,keep_offset);
saveas(gcf, strcat(imagename, 'GCR_AIC_flexible'),'epsc')

figure;
make_plot_GCR_color_polygon(X,fourier_coef,fourier_coef_b_subset,[1,2],candidate_bic,1000,colors);
gca = image_rescaling(gca,dataname,uncertain_plot_include_x,uncertain_plot_include_y,keep_offset);
saveas(gcf, strcat(imagename, 'GCR_BIC_flexible'),'epsc')

%% plot pdf 
pdf_figuresize_simple = [50 50 450 320];
legdvalue = {'bootstrap','obeserved'};
fontsize = 20;
legdsize = 15;
ticks_num=6;
legdpos = 'eastoutside';
%% plot flux (in bootstrap)

lwidth=2;
drop = false;
binnum = 10;
relfreq = true;
linetypes = {'-','-.','--'}; linecolors = {'k','k','k'};
num_in_poly = sum(inpolygon(X(:,1),X(:,2),contourx,contoury));
area_poly = polyarea(contourx,contoury);
num_in_adjust_poly = get_adjust_flux(num_in_poly,area_poly,n);

figure
flux_btp = {};
flux_btp{1} = get_adjust_flux(num_in_poly_set,area_poly_set,n);

make_freq_stack(flux_btp,[num_in_adjust_poly], ...
    {}, legdpos,linetypes,linecolors,binnum,lwidth,fontsize,drop,legdsize,relfreq,21)
xlabel('Adjust Flux')
%ylabel('f(Adjust Flux)')
set(gcf, 'Position', pdf_figuresize_simple); %
saveas(gcf, strcat(imagename, 'bootstrap_adjustflux_simple'),'epsc')
quantiles = quantile(flux_btp{1} , [0.16, 0.84]);
mean_val = mean(flux_btp{1});

% Compute the mode
mode_val = mode(flux_btp{1});
% Display the results
fprintf('adjust flux 16th quantile: %.2f\n', quantiles(1));
fprintf('adjust flux 84th quantile: %.2f\n', quantiles(2));
fprintf('adjust flux Mean: %.2f\n', mean_val);
fprintf('adjust flux Mode: %.2f\n', mode_val);
fprintf('adjust flux: %.2f\n', num_in_adjust_poly);

%% plot area (in bootstrap)

figure
area_btp = {};
area_btp{1} = area_poly_set*convertrate;
area_binwith = 700;
area_binnum = round((max(area_btp{1})-min(area_btp{1}))/area_binwith);
area_poly_convert = area_poly*convertrate;
make_freq_stack(area_btp,[area_poly_convert], ...
    {}, legdpos,linetypes,linecolors,area_binnum,lwidth,fontsize,drop,legdsize,relfreq,21,[0,9000])
xlabel('Area [arcsec^2]',"linewidth",2)
%ylim( [0,0.22] )
set(gcf, 'Position', pdf_figuresize_simple); %
saveas(gcf, strcat(imagename, 'bootstrap_area_simple'),'epsc')
quantiles = quantile(area_btp{1}  , [0.16, 0.84]);
mean_val = mean(area_btp{1} );

% Compute the mode
mode_val = mode(area_btp{1} );
% Display the results
fprintf('area 16th quantile: %.2f\n', quantiles(1));
fprintf('area 84th quantile: %.2f\n', quantiles(2));
fprintf('area Mean: %.2f\n', mean_val);
fprintf('area Mode: %.2f\n', mode_val);
fprintf('area: %.2f\n', area_poly_convert);



%% plot the centroid bootstrap
figure
centx = real(fourier_coef(1));centy = imag(fourier_coef(1)); 
centx_all = real(fourier_coef_b_subset(1,:));centy_all = imag(fourier_coef_b_subset(1,:));


make_centroid_plot_with_hist([centx,centy],[centx_all;centy_all]','boxplot')
gca = image_rescaling(gca,strcat(dataname,'centroid'));
saveas(gcf, strcat(imagename, 'centroid_scatter_box'),'epsc')

figure
make_centroid_plot_with_hist([centx,centy],[centx_all;centy_all]','histogram')
gca = image_rescaling(gca,strcat(dataname,'centroid'));
saveas(gcf, strcat(imagename, 'centroid_scatter_hist'),'epsc')

%% get the centroid and std 

calc_arcra = @(x) calc_arcsec(x,4800.0,23707.5,26107.5,0.05,5.38516,1);
calc_arcdec = @(x) calc_arcsec(x,4800.0,25646.5,28046.5,0.05,-6.92139,2);
x_std_val = std(calc_arcra(real(fourier_coef_b_subset(1,:))));
y_std_val = std(calc_arcdec(imag(fourier_coef_b_subset(1,:))));
x_val = calc_arcra(real(fourier_coef(1)));
y_val = calc_arcdec(imag(fourier_coef(1)));
fprintf('X centroid: %.2f\n', x_val);
fprintf('Y centroid: %.2f\n', y_val);
fprintf('X centroid Standard Deviation: %.2f\n', x_std_val);
fprintf('Y centroid Standard Deviation: %.2f\n', y_std_val);
figure
scatter(calc_arcra(real(fourier_coef_b_subset(1,:))),calc_arcdec(imag(fourier_coef_b_subset(1,:))),'.')
hold on
scatter(calc_arcra(real(fourier_coef(1))),calc_arcdec(imag(fourier_coef(1))),'filled')
set(gca, 'xdir', 'reverse');
axis equal
xlim([-30 30])
ylim([-30 30])

%% plot the eccentricity 
%rescale to real data scale

fourier_coef_b_subset_rescale = rescale_coef(fourier_coef_b_subset,0,0,convertrate_x,convertrate_y);
long_all_axis = abs(abs(fourier_coef_b_subset_rescale(2,:))+abs(fourier_coef_b_subset_rescale(end,:))); 
short_all_axis = abs(abs(fourier_coef_b_subset_rescale(2,:))-abs(fourier_coef_b_subset_rescale(end,:))); 
fourier_coef_rescale = rescale_coef(fourier_coef,0,0,convertrate_x,convertrate_y);
long_axis = abs(abs(fourier_coef_rescale(2))+abs(fourier_coef_rescale(end))); 
short_axis = abs(abs(fourier_coef_rescale(2))-abs(fourier_coef_rescale(end))); 
figure
make_eccentricity_plot_with_hist([long_axis,short_axis],[long_all_axis;short_all_axis]',minor_axis_limit,major_axis_limit, dataname,'boxplot')
saveas(gcf, strcat(imagename, 'eccentricity_scatter_box'),'epsc')
figure
make_eccentricity_plot_with_hist([long_axis,short_axis],[long_all_axis;short_all_axis]',minor_axis_limit,major_axis_limit, dataname,'histogram')
saveas(gcf, strcat(imagename, 'eccentricity_scatter_hist'),'epsc')


%% get eccentricity and std
[invx,invy] = iFD  (fourier_coef);
convertx = calc_arcra(invx);
converty = calc_arcdec(invy);
fourier_coef_rescale2 = FD(convertx,converty);
% x_val = abs(fourier_coef_rescale2(2))+abs(fourier_coef_rescale2(end));
% y_val = abs(abs(fourier_coef_rescale2(2))-abs(fourier_coef_rescale2(end)));
x_val = long_axis;
y_val = short_axis;
x_std_val = std(long_all_axis);
y_std_val = std(short_all_axis);


fprintf('major axis: %.2f\n', x_val);
fprintf('minor axis: %.2f\n', y_val);
fprintf('major axis Standard Deviation: %.2f\n', x_std_val);
fprintf('minor axis Standard Deviation: %.2f\n', y_std_val);

%% plot 10 illustration of bootstrap results
figure 
contourx_all_subset = contourx_all(~invalid_index);
contoury_all_subset = contoury_all(~invalid_index);

hold on
for ii = 1:10
    contourx_b = contourx_all_subset{ii};
    contoury_b = contoury_all_subset{ii};
    plot(contourx_b,contoury_b,'Color',GRAY,LineWidth=2)
end
box on 
axis equal
gca = image_rescaling(gca,dataname,uncertain_plot_include_x,false,keep_offset);
plot(contourx,contoury,'r',LineWidth=2)
saveas(gcf, strcat(imagename, 'illustration_bootstrap_10'),'epsc')

