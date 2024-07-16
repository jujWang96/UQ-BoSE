parpath = '~/src/gsrg/';
plotpath = 'Astro_sim/real_data/plots/';
resultpath = 'Astro_sim/real_data/results/';
addpath(genpath(strcat(parpath,'Astro_sim')))
addpath(genpath(strcat(parpath,'G-SRG')))
addpath(genpath(strcat(parpath,'plot_util')))
addpath(genpath(strcat(parpath,'contour_util')))
addpath(genpath(strcat(parpath,'FD_util')))
addpath(genpath(strcat(parpath,'boot_util')))
addpath(genpath(strcat(parpath,'gcr_util')))
addpath(genpath(strcat(parpath,'Arp299XMM')))

dataname = "Arp299XMM";

load(strcat(parpath,resultpath,dataname,"10000rep_vor_UQ_analysis.mat"))

imagename = strcat(parpath, plotpath, dataname, "_10000rep_vor_UQ_");

%% plot parameters
line_size = 2;

convertrate_x = 35.2;
convertrate_y = 35.2;
convertrate = convertrate_x*convertrate_y*100;
major_axis_limit = 50;
minor_axis_limit = 10;
keep_offset = false;
uncertain_plot_include_x = false;
uncertain_plot_include_y = false;

[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
[V, R] = voronoiDiagram(DT);


%% plot the histogram of voronoi area 
figure
histogram(cell_area)

%% scatter plot 
figure;
scatter(X(:,1),X(:,2),1,'k.')
gca = image_rescaling(gca,dataname,true,true,keep_offset);

% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
%[seeds, num] = get_seeds_sim_kde(cx, cy,invalid, adj_mat, threshold);
[seeds, num] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat,cell_area,best_P, best_threshold);
disp(strcat("Num of seeds is ",num2str(num)))
disp(strcat("Value of FD (BIC) is ",num2str(candidate_bic)))
disp(strcat("Number of segments is ",num2str(length(raw_contour.regionId))))



%% plot contours 
figure;

target = [contourx,contoury];
make_plot_contour(X,raw_contour,1,line_size,target);
[invx,invy] = iFD(fourier_coef,candidate_bic);
hold on
plot(invx,invy,'r',LineWidth=line_size)
gca = image_rescaling(gca,dataname,true,false,keep_offset);
saveas(gcf, strcat(imagename, 'disconnect_boundary'),'epsc')


%% plot the flexible model selection results 
figure ;
make_plot_modelselect_flexible(fourier_coef,candidate_aic,pgon_confine,line_size);
gca = image_rescaling(gca,dataname,uncertain_plot_include_x,uncertain_plot_include_y,keep_offset);

figure;
make_plot_modelselect_flexible(fourier_coef,candidate_bic,pgon_confine,line_size);
gca = image_rescaling(gca,dataname,uncertain_plot_include_x,uncertain_plot_include_y,keep_offset);
%saveas(gcf, strcat(imagename, 'flexible_model_selection_BIC'),'epsc')
