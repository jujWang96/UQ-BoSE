

load('model_select_result_2300chandra.mat')
load('ngc2300_box_058kev_evt2.mat')
addpath(genpath('~/Desktop/gsrg/boot_util'))
target_idx = 4;
%P = 14; threshold = 15;
X = double(unique(X,'rows'));
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(double(X), [0 1], [0 1], ones(size(X, 1), 1));

min_pair_num_AIC = min_pair_num{1};
min_select_AIC = min_pair_num_AIC(target_idx)*2+1;

min_pair_num_BIC = min_pair_num{2};
min_select_BIC = min_pair_num_BIC(target_idx)*2+1;
smooth_polyset_AIC = smooth_polyset{1};
smooth_polyset_BIC = smooth_polyset{2};
B = 1;

%initialize for AIC
fourier_coef_all = zeros(B,300);
contourx_all = cell(1,B);
contoury_all = cell(1,B);
area_poly_all = zeros(1,B);
area_aic_all = zeros(1,B);
num_B_poly_all = zeros(1,B); 
num_B_aic_all = zeros(1,B);
intensity_B_poly_all = zeros(1,B);
intensity_B_aic_all = zeros(1,B);
num_O_poly_all = zeros(1,B); 
num_O_aic_all = zeros(1,B); 
intensity_O_poly_all = zeros(1,B);
intensity_O_aic_all = zeros(1,B);
centx_poly_all = zeros(1,B);
centy_poly_all = zeros(1,B);
centx_aic_all = zeros(1,B); 
centy_aic_all = zeros(1,B);
%initialize for BIC
fourier_coef_all_BIC = zeros(B,300);
contourx_all = cell(1,B);
contoury_all = cell(1,B);
area_poly_all = zeros(1,B);
area_aic_all = zeros(1,B);
num_B_poly_all = zeros(1,B); 
num_B_aic_all = zeros(1,B);
intensity_B_poly_all = zeros(1,B);
intensity_B_aic_all = zeros(1,B);
num_O_poly_all = zeros(1,B); 
num_O_aic_all = zeros(1,B); 
intensity_O_poly_all = zeros(1,B);
intensity_O_aic_all = zeros(1,B);
centx_poly_all = zeros(1,B);
centy_poly_all = zeros(1,B);
centx_aic_all = zeros(1,B); 
centy_aic_all = zeros(1,B);
parfor seed = 1:B
    [bootstrap_X_set,new_seeds_set] = sim_fix_data_vor_withseed(X,cell_area, n, DT,2,1,seed,seeds_all);
    [fourier_coef_all(seed,:), contourx_all{seed},contoury_all{seed}, area_poly_all(seed),area_aic_all(seed), area_bic_all(seed), num_B_poly_all(seed), num_B_aic_all(seed),num_B_bic_all(seed), intensity_B_poly_all(seed),...
    intensity_B_aic_all(seed), intensity_B_bic_all(seed), num_O_poly_all(seed),num_O_aic_all(seed), num_O_bic_all(seed), intensity_O_poly_all(seed), intensity_O_aic_all(seed),intensity_O_bic_all(seed),...
    centx_poly_all(seed), centy_poly_all(seed), centx_aic_all(seed), centy_aic_all(seed),centx_bic_all(seed), centy_bic_all(seed)] = sim_fit_random_real_data_withseed(bootstrap_X_set{1}, seed, false, P, threshold, ...
    rand_num, rep_itr, X, min_select_AIC,min_select_BIC, penalty,M, new_seeds_set{1}); 

end

save('bootstrap_result_2300chandra_vor_withseed.mat')