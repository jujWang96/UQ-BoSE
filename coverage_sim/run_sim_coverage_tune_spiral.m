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
GRAY = [0.6 0.6 0.6];

Beta = 0.25:0.25:3;

p = 1-exp(-Beta.^2/2);
%% srgong parameter
gridsize = 0.1;
seedsize = 5;
rep_itr = 1000;

grid_only = true;
force_merge = true;
show_plot = false;
merge_method = "random";
r = 1000;


%% data generate parameter
shapename = 'spiral';
true_FD_num = 11;
V(:,1) = V_spiral(:,1);
V(:,2) = V_spiral(:,2);
factor = 30; %factors = [30,60];
sample_factor = 1; %sample_factors = [1,2];
lambda =1000;
base_num_in_circle = 10 ;

%% simulation parameters
B = 200;
numset = 500;

%% model selection and bootstrap parameter 
selection_method = "true_FD_num"%"BIC";
pair_search = true;

boot_method = "true_FD_num"%"BIC";
K = 300;


origin_coef_set = zeros(K,numset);

covers = zeros(numset, length(p));
GCRs = {};        
totalN = (factor*base_num_in_circle+lambda)*sample_factor;

parfor rep = 1:numset
    try
        rng(rep,'twister')
        objN = binornd(totalN,sample_factor*factor * base_num_in_circle/totalN);
        %rep
        [GCRs{rep}, covers(rep,:),origin_coef_set(:,rep)]= sim_coverage_tune_get_prob(V,totalN, objN, rep,gridsize,seedsize, grid_only, merge_method,force_merge, rep_itr, pair_search, selection_method, boot_method,B, p, true_FD_num);
    catch
        
         error(['Error in File ' num2str(rep)])

    end
end

invalid_index = all(origin_coef_set==0);
%covers(invalid_index,:) = [];
covers = covers(~invalid_index,:);
origin_coef_set_truncate = origin_coef_set([end-div2(true_FD_num-1)+1:end,1:div2(true_FD_num)+1],:);

origin_coef_set_truncate = origin_coef_set_truncate( :,~invalid_index); 

para = transpose([real(origin_coef_set_truncate);imag(origin_coef_set_truncate)]);
para(:,true_FD_num+div2(true_FD_num-1)+2)=[];
C = cov(para);
p_LCs = zeros(1,length(Beta));
p_Gs = zeros(1,length(Beta));
for ii = 1:length(Beta)
    p_Gs(ii) =  chi2cdf(Beta(ii)^2,2*true_FD_num-1);
end
p_LCs = max(0,real(calc_lc_bd_single_continuous(true_FD_num,C,Beta)));

p_sims = mean(covers,1);

save(strcat(parpath,resultpath,"coverage_prob_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor)  , "_gridsize", num2str(gridsize), ...
            "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search),...
            "_selection_method",selection_method, "_boot_method" ,boot_method,'.mat'))
