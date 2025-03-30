clear 
close all
parpath = '~/src/gsrg/';
resultpath = 'coverage_sim/results/';
plotpath = 'coverage_sim/plots/';
addpath(genpath(strcat(parpath,'coverage_sim')))

addpath(genpath(strcat(parpath,'Astro_sim')))

addpath(genpath(strcat(parpath,'G-SRG')))
addpath(genpath(strcat(parpath,'plot_util')))
addpath(genpath(strcat(parpath,'contour_util')))
addpath(genpath(strcat(parpath,'FD_util')))
addpath(genpath(strcat(parpath,'boot_util')))
addpath(genpath(strcat(parpath,'gcr_util')))
addpath(genpath(strcat(parpath,'Astro_sim/sim')))
addpath(genpath(strcat(parpath,'Astro_sim/sim_util')))

shapenames = ["ellipse","bridge", "spiral"] ;
true_FD_nums =[3,7,11];


r = 1000;
GRAY = [0.6 0.6 0.6];


%% data generate parameter

factor = 30; %factors = [30,60];
sample_factor = 1; %sample_factors = [1,2];
lambda =1000;
base_num_in_circle = 10 ;

%% srgong parameter
gridsize = 0.1;
seedsize = 5;
grid_only = true;
force_merge = true;
pair_search = true;
show_plot = true;
merge_method = "random";
rep_itr = 1000;

B = 200;
K = 300;
numset = 300;
subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.06], [0.1 0.02], [0.08 0.02]);
candidate_bics_filter = {};
candidate_aics_filter = {};
for  ii = 1: length(shapenames)
    shapename = shapenames(ii);
    true_FD_num = true_FD_nums(ii);
    
    load(strcat(parpath,resultpath, "model_selection_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor)  , "_gridsize", num2str(gridsize),...
                "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search),'.mat'))

%                     load(strcat("model_selection_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor)  ,...
%                     "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search),'.mat'))
        candidate_bics_filter{ii} = candidate_bics(candidate_bics ~= 0);
        candidate_aics_filter{ii} = candidate_aics(candidate_aics ~= 0);
%         figure
%         histogram_plot_group(candidate_aics_filter, true_FD_num, factor, sample_factor)
    
end
figure
histogram_plot_group_model_selection(candidate_bics_filter, true_FD_nums,factor, sample_factor,true)
set(gcf,'position',[0,0,600,180])
saveas(gcf,strcat(parpath,plotpath, "BIC","_sample_factor", num2str(sample_factor),"_factor", num2str(factor)  , "_gridsize", num2str(gridsize),...
                "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search)),'epsc')
figure
histogram_plot_group_model_selection(candidate_aics_filter, true_FD_nums,factor, sample_factor,false)
set(gcf,'position',[0,0,600,180])
saveas(gcf,strcat(parpath,plotpath, "AIC","_sample_factor", num2str(sample_factor),"_factor", num2str(factor)  , "_gridsize", num2str(gridsize),...
                "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search)),'epsc')
