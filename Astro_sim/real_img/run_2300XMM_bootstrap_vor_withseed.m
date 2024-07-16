

load('model_select_result_2300XMM.mat')
load('ngc2300_MOS1_evt_0.5-8.0keV_scaled.mat')
X = double(unique(X,'rows'));
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));

%P = 14; threshold=20;
target_idx = 3;
min_pair_num_AIC = min_pair_num{1};
min_select_AIC = min_pair_num_AIC(target_idx)*2+1;

min_pair_num_BIC = min_pair_num{2};
min_select_BIC = min_pair_num_BIC(target_idx)*2+1;
smooth_polyset_AIC = smooth_polyset{1};
smooth_polyset_BIC = smooth_polyset{2};
B = 200;

%initialize for AIC
fourier_coef_all_AIC = zeros(B,300);
contourx_all_AIC = cell(1,B);
contoury_all_AIC = cell(1,B);
area_poly_all_AIC = zeros(1,B);
area_fd_all_AIC = zeros(1,B);
num_B_poly_all_AIC = zeros(1,B); 
num_B_fd_all_AIC = zeros(1,B);
intensity_B_poly_all_AIC = zeros(1,B);
intensity_B_fd_all_AIC = zeros(1,B);
num_O_poly_all_AIC = zeros(1,B); 
num_O_fd_all_AIC = zeros(1,B); 
intensity_O_poly_all_AIC = zeros(1,B);
intensity_O_fd_all_AIC = zeros(1,B);
centx_poly_all_AIC = zeros(1,B);
centy_poly_all_AIC = zeros(1,B);
centx_fd_all_AIC = zeros(1,B); 
centy_fd_all_AIC = zeros(1,B);
%initialize for BIC
fourier_coef_all_BIC = zeros(B,300);
contourx_all_BIC = cell(1,B);
contoury_all_BIC = cell(1,B);
area_poly_all_BIC = zeros(1,B);
area_fd_all_BIC = zeros(1,B);
num_B_poly_all_BIC = zeros(1,B); 
num_B_fd_all_BIC = zeros(1,B);
intensity_B_poly_all_BIC = zeros(1,B);
intensity_B_fd_all_BIC = zeros(1,B);
num_O_poly_all_BIC = zeros(1,B); 
num_O_fd_all_BIC = zeros(1,B); 
intensity_O_poly_all_BIC = zeros(1,B);
intensity_O_fd_all_BIC = zeros(1,B);
centx_poly_all_BIC = zeros(1,B);
centy_poly_all_BIC = zeros(1,B);
centx_fd_all_BIC = zeros(1,B); 
centy_fd_all_BIC = zeros(1,B);
parfor seed = 1:B
    [bootstrap_X_set,new_seeds_set] = sim_fix_data_vor_withseed(X,cell_area, n, DT,2,1,seed,seeds_all);
    [fourier_coef_all(seed,:), contourx_all{seed},contoury_all{seed}, area_poly_all(seed),area_aic_all(seed), area_bic_all(seed), num_B_poly_all(seed), num_B_aic_all(seed),num_B_bic_all(seed), intensity_B_poly_all(seed),...
    intensity_B_aic_all(seed), intensity_B_bic_all(seed), num_O_poly_all(seed),num_O_aic_all(seed), num_O_bic_all(seed), intensity_O_poly_all(seed), intensity_O_aic_all(seed),intensity_O_bic_all(seed),...
    centx_poly_all(seed), centy_poly_all(seed), centx_aic_all(seed), centy_aic_all(seed),centx_bic_all(seed), centy_bic_all(seed)] = sim_fit_random_real_data_withseed(bootstrap_X_set{1}, seed, false, P, threshold, ...
    rand_num, rep_itr, X, min_select_AIC,min_select_BIC, penalty,M, new_seeds_set{1}); 

end

save('bootstrap_result_2300XMM_vor_withseed.mat')