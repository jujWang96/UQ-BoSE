addpath(genpath('/home/minjay/G-SRG'))
addpath(genpath('/home/minjay/Astro_sim'))

load('real_full_result_2018_11_1_23_4_37.mat')

n_boot = 100;
impute_log_int = zeros(length(cx), n_boot);
n_region = zeros(n_boot, 1);
parfor i = 1:n_boot
    [impute_log_int(:, i), n_region(i)] = real_full_jitter(i);
end

filename = 'real_full_result_boot.mat';
save(filename, 'impute_log_int', 'n_region')