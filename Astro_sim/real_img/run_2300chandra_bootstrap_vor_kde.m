

load('model_select_result_2300chandra.mat')
load('ngc2300_box_058kev_evt2.mat')
addpath(genpath('~/Desktop/gsrg/boot_util'))
target_idx = 4;
%P = 14; threshold = 15;
X = double(unique(X,'rows'));

[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));


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
contourx_all = cell(1,B);
contoury_all = cell(1,B);
area_poly_all = zeros(1,B);
area_fd_all = zeros(1,B);
num_B_poly_all = zeros(1,B); 
num_B_fd_all = zeros(1,B);
intensity_B_poly_all = zeros(1,B);
intensity_B_fd_all = zeros(1,B);
num_O_poly_all = zeros(1,B); 
num_O_fd_all = zeros(1,B); 
intensity_O_poly_all = zeros(1,B);
intensity_O_fd_all = zeros(1,B);
centx_poly_all = zeros(1,B);
centy_poly_all = zeros(1,B);
centx_fd_all = zeros(1,B); 
centy_fd_all = zeros(1,B);
gridx1 = 0.005:.005:0.995;
gridx2 = 0.005:.005:0.995;

[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
[cdf_values, cdf_x] = mvksdensity(X, xi,'Kernel', 'epanechnikov','Function','cdf','support',[-0.001,-0.001;1.001,1.001]);
n_samples = length(X);
parfor seed = 1:B
    bootstrap_X_set_AIC = generate_samples_from_kde(X,cdf_values, cdf_x,gridx1,gridx2, n_samples,seed);
    [fourier_coef_all_AIC(seed,:), contourx_all_AIC{seed},contoury_all_AIC{seed}, area_poly_all_AIC(seed),area_fd_all_AIC(seed), num_B_poly_all_AIC(seed), num_B_fd_all_AIC(seed), intensity_B_poly_all_AIC(seed),...
    intensity_B_fd_all_AIC(seed), num_O_poly_all_AIC(seed),num_O_fd_all_AIC(seed), intensity_O_poly_all_AIC(seed), intensity_O_fd_all_AIC(seed),...
    centx_poly_all_AIC(seed), centy_poly_all_AIC(seed), centx_fd_all_AIC(seed), centy_fd_all_AIC(seed)] = sim_fit_random_real_data(bootstrap_X_set_AIC, seed, false, P, threshold, ...
    rand_num, rep_itr, X, min_select_AIC,penalty,M);
    bootstrap_X_set_BIC = generate_samples_from_kde(X,cdf_values, cdf_x,gridx1,gridx2, n_samples,seed);

    [fourier_coef_all_BIC(seed,:), contourx_all_BIC{seed},contoury_all_BIC{seed}, area_poly_all_BIC(seed),area_fd_all_BIC(seed), num_B_poly_all_BIC(seed), num_B_fd_all_BIC(seed), intensity_B_poly_all_BIC(seed),...
    intensity_B_fd_all_BIC(seed), num_O_poly_all_BIC(seed),num_O_fd_all_BIC(seed), intensity_O_poly_all_BIC(seed), intensity_O_fd_all_BIC(seed),...
    centx_poly_all_BIC(seed), centy_poly_all_BIC(seed), centx_fd_all_BIC(seed), centy_fd_all_BIC(seed)] = sim_fit_random_real_data(bootstrap_X_set_BIC, seed, false, P, threshold, ...
    rand_num, rep_itr, X, min_select_BIC,penalty,M);

end

save('bootstrap_result_2300chandra_kde.mat')
