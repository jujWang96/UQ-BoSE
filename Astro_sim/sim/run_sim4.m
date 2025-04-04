addpath(genpath('~/UQ-BoSE/Astro_sim'))
addpath(genpath('~/UQ-BoSE/G-SRG'))
addpath(genpath('~/UQ-BoSE/plot_util'))

loc = [0.5 0.5; 0.4 0.4; 0.4 0.6; 0.6 0.4; 0.6 0.6];
radius = [0.25 0.025*ones(1, 4)];
base_num_in_circle = [10 ones(1, 4)];
factors = [10 20 30];
lambda = 1000;
sample_factors = [0.5 1 2];
T = 500;
metrics = zeros(length(factors) * length(sample_factors), T);
n_region = zeros(length(factors) * length(sample_factors), T);
merge_time = zeros(length(factors) * length(sample_factors), T);
gridszie = 0.2;
seedsize = 5;
for seed = 1:T
        [metrics(:, seed), n_region(:, seed), merge_time(:,seed)] = sim4(loc, radius, base_num_in_circle, factors, lambda, sample_factors, seed,gridszie, seedsize);
end

save(strcat('sim4_result_seedsize5_gs',num2str(gridszie),'.mat'))
   
metrics = zeros(length(factors) * length(sample_factors), T);
n_region = zeros(length(factors) * length(sample_factors), T);
merge_time = zeros(length(factors) * length(sample_factors), T);

parfor seed = 1:T
        [metrics(:, seed), n_region(:, seed), merge_time(:,seed)] = sim4(loc, radius, base_num_in_circle, factors, lambda, sample_factors, seed,gridszie);
end

save(strcat('sim4_result_gs',num2str(gridszie),'.mat'))
   
