close all
clear
load('model_select_result_2300chandra.mat')
load('ngc2300_box_058kev_evt2.mat')
addpath(genpath('~/Desktop/gsrg/boot_util'))
cd '/Users/joycewang/Desktop/gsrg/NGC2300'
target_idx = 4;
X = double(unique(X,'rows'));
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));

%P = 14; threshold = 15;
min_pair_num_AIC = min_pair_num{1};
min_select_AIC = min_pair_num_AIC(target_idx)*2+1;

min_pair_num_BIC = min_pair_num{2};
min_select_BIC = min_pair_num_BIC(target_idx)*2+1;
smooth_polyset_AIC = smooth_polyset{1};
smooth_polyset_BIC = smooth_polyset{2};
B = 1;


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
origin_invx = target_x{target_idx};
origin_invy = target_y{target_idx};
threshold = 15
boot_method = 2

% %plot the histogram of 10 bootstrapped data 
% boot_X_mult= [];
% r = 100;
% for seed = 1:r
%     [bootstrap_X_set,~] = sim_fix_data_vor_fast(X,cell_area, n, DT,3,1,seed);
%     boot_X_mult = [boot_X_mult;bootstrap_X_set{1}];
% end
% sample_r = randsample(r*length(X),length(X));
% boot_X_mult_sub = boot_X_mult(sample_r,:);
% boot_X_mult_sub = double(unique(boot_X_mult_sub,'rows'));
% [~, ~, ~, ~, ~, ~, cell_area_boot] = init_comp(boot_X_mult_sub, [0 1], [0 1], ones(size(boot_X_mult_sub, 1), 1));
% figure
% histogram(cell_area, "FaceAlpha", 0.5)
% hold on
% histogram(cell_area_boot,"FaceAlpha", 0.5)
% saveas(gcf, 'chandra2300_bootstrap3.png')

% boot_X_mult= [];
% r = 100;
% for seed = 1:r
%     [bootstrap_X_set,~] = sim_fix_data_vor_fast(X,cell_area, n, DT,1,1,seed);
%     boot_X_mult = [boot_X_mult;bootstrap_X_set{1}];
% end
% sample_r = randsample(r*length(X),length(X));
% boot_X_mult_sub = boot_X_mult(sample_r,:);
% boot_X_mult_sub = double(unique(boot_X_mult_sub,'rows'));
% [~, ~, ~, ~, ~, ~, cell_area_boot] = init_comp(boot_X_mult_sub, [0 1], [0 1], ones(size(boot_X_mult_sub, 1), 1));
% figure
% histogram(cell_area, "FaceAlpha", 0.5)
% hold on
% histogram(cell_area_boot,"FaceAlpha", 0.5)
% saveas(gcf, 'chandra2300_bootstrap1.png')

%P=14
%rep_itr=5000
%M=200
for seed = 1:10
    [bootstrap_X_set_AIC,Regions] = sim_fix_data_vor_fast(X,cell_area, n, DT,boot_method,1,seed);
    sim_fit_2300chandra_bootstrap(bootstrap_X_set_AIC{1}, seed, false, P, threshold, ...
    rand_num, rep_itr, X, min_select_AIC,penalty,M,origin_invx,origin_invy,boot_method);
    
end
