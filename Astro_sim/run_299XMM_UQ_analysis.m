close all;
clear;
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

load('Arp299_MOS1_evt_0.5-8.0keV_scaled.mat')
dataname = "Arp299XMM";
X = double(unique(X,'rows'));
n = length(X);
K =300;
B = 200;
r = 1024;
colors = [[0.6,0.6,0.8];[0.8,0.8,1]];

GRAY = [0.6,0.6,0.6];
%% seeds parameter
Ps=[5,8,11,14,17,20];
thresholds = [5,10,15,20];
%% srgong parameter to tune
M = 150;
rand_num = 5;
rep_itr = 10000;
penalty = 6;
colors = [[0.6,0.6,0.8];[0.8,0.8,1]];

%% perform seeds parameter tuning for P and threshold
BIC_tune = zeros(1,length(Ps)*length(thresholds));
BIC_tune_2seg = zeros(1,length(Ps)*length(thresholds));
best_P = 0;
best_threshold = 0;
best_BIC= inf;
idx = 0;
for P = Ps
    for threshold = thresholds
        idx = idx+1;
        [contourx,contoury, fourier_coef, raw_contour,selected,min_BIC, min_BIC_2seg] = sim_fit_random_stratified(X, idx, false, P, threshold, rand_num, rep_itr,penalty, M);
        BIC_tune(idx) = min_BIC;
        BIC_tune_2seg(idx) = min_BIC_2seg;
        if min_BIC_2seg < best_BIC
            best_P = P;
            best_threshold = threshold;
            best_BIC = min_BIC_2seg;
        end
    end
end
disp(strcat("The best seed parameter is", "P=",num2str(best_P), ", threshold=",num2str(best_threshold)))
[~,param_idx] = min(BIC_tune_2seg);
P = Ps(ceil(param_idx/length(thresholds)));
threshold = thresholds((mod(param_idx-1,length(thresholds)))+1);

%% fit model with tuned parameter and perform model selection 
[contourx,contoury, fourier_coef, raw_contour,selected, min_BIC, min_BIC_2seg] = sim_fit_random_stratified(X, idx, false, best_P, best_threshold, rand_num, rep_itr,penalty, M);

plot_results = true;
candidates = 3:2:K ;
[candidate_aic,candidate_bic,pgon_confine] = model_selection_flexible_robust(X,contourx, contoury, raw_contour,candidates,plot_results);

%% generate B nonparametric bootstrap data and derive uncertainty statistics
fourier_coef_b_set = zeros(K,B);
area_poly_set = zeros(1,B);
num_in_poly_set = zeros(1,B);
intensity_poly = zeros(1,B);
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));

parfor seed = 1:B
    bootstrap_X_set_vor = sim_fix_data_vor_fast_vk(X,cell_area, n, DT,1,seed);
    bootstrap_X = bootstrap_X_set_vor{1};
    [contourx_b,contoury_b, fourier_coef_b, raw_contour_b,selected_nonempty_b, min_BIC_b, min_BIC_2seg_b] = sim_fit_random_stratified(bootstrap_X, seed, false, best_P, best_threshold, rand_num, rep_itr,penalty, M);
    fourier_coef_b_set(:,seed) = fourier_coef_b;
    contourx_all{seed} = contourx_b;
    contoury_all{seed} = contoury_b;
    area_poly_set(seed) = polyarea(contourx_b,contoury_b);
    num_in_poly_set(seed) = sum(inpolygon(bootstrap_X(:,1),bootstrap_X(:,2),contourx_b,contoury_b));

end
delete(gcp)
invalid_index = all(fourier_coef_b_set==0);
fourier_coef_b_subset = fourier_coef_b_set( :,~invalid_index);
make_plot_GCR_color_polygon(X,fourier_coef,fourier_coef_b_subset,[1,2],candidate_bic,1000);
gca = image_rescaling(gca,dataname,true,true,true);

save(strcat(parpath,resultpath,dataname,"10000rep_vor_UQ_analysis.mat"))

