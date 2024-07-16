close all
cd /Users/joycewang/Desktop/gsrg/Arp299XMM/

load('tune_param_Arp299.mat')
load('model_select_result_299XMM.mat')

load('Arp299_MOS1_evt_0.5-8.0keV_scaled.mat')
X = double(unique(X,'rows'));

min_pair_num_AIC = min_pair_num{1};

addpath(genpath('~/Desktop/gsrg/boot_util'))


[fourier_coef, contourx,contoury, area_poly,area_fd, num_B_poly, num_B_fd, intensity_B_poly,...
    intensity_B_fd, num_O_poly,num_O_fd, intensity_O_poly, intensity_O_fd,...
    centx_poly, centy_poly, centx_fd, centy_fd] = sim_fit_random_real_data(X, 0, false, P, threshold, rand_num, rep_itr, X, min_pair_num_AIC(target_idx)*2+1,penalty, M);
