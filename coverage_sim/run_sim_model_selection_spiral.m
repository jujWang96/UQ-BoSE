clear 
close all
parpath = '~/UQ-BoSE/';
resultpath = 'coverage_sim/results/';
addpath(genpath(strcat(parpath,'Astro_sim')))
addpath(genpath(strcat(parpath,'G-SRG')))
addpath(genpath(strcat(parpath,'plot_util')))
addpath(genpath(strcat(parpath,'contour_util')))
addpath(genpath(strcat(parpath,'FD_util')))
addpath(genpath(strcat(parpath,'boot_util')))
addpath(genpath(strcat(parpath,'gcr_util')))
addpath(genpath(strcat(parpath,'Astro_sim/sim')))
addpath(genpath(strcat(parpath,'Astro_sim/sim_util')))

load('shape_coordinates.mat')

true_FD_num = 11;
r = 1000;
GRAY = [0.6 0.6 0.6];

Beta = 0.25:0.25:4;

p = 1-exp(-Beta.^2/2);

%% data generate parameter
shapename = 'spiral';
V(:,1) = V_spiral(:,1);
V(:,2) = V_spiral(:,2);
factors = [30,40];
sample_factors = [1,2,3];
lambda =1000;
base_num_in_circle = 10 ;

%% srgong parameter
merge_method = "random";

gridsize = 0.1;
seedsize = 5;
rep_itr = 1000;

grid_only = true;
force_merge = true;
pair_search = false;
robust = true;
show_plot = true;

B = 200;
K = 300;
numset = 300;

for factor = factors
    for sample_factor = sample_factors
        totalN = (factor*base_num_in_circle+lambda)*sample_factor;
        candidate_aics = zeros(1,numset);
        candidate_bics = zeros(1,numset);
        parfor rep = 1:numset
            try
	            rng(rep,'twister')
                objN = binornd(totalN,sample_factor*factor * base_num_in_circle/totalN);
            
                [candidate_aics(rep),candidate_bics(rep)] = sim_coverage_tune_model_selection(V,totalN, objN, rep,gridsize,seedsize, grid_only, merge_method,force_merge, rep_itr, pair_search, robust);
                
                %origin_coef_set(:,rep) = origin_coef_inv([end-div2(M-1)+1:end,1:div2(M)+1]);
    	    catch
	   	        error(['Error in File ' num2str(rep)])
	        end    
        end
        save(strcat(parpath,resultpath, "model_selection_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor)  , "_gridsize", num2str(gridsize), ...
            "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search),'.mat'))
delete(gcp)

    end
end
