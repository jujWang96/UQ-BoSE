addpath(genpath('/home/minjay/G-SRG'))
addpath(genpath('/home/minjay/Astro_sim'))

loc = [0.5 0.5; 0.4 0.4; 0.4 0.6; 0.6 0.4; 0.6 0.6];
radius = [0.25 0.025*ones(1, 4)];
base_num_in_circle = [10 ones(1, 4)];
factors = [10 20 30];
lambda = 1000;
sample_factors = [0.5 1 2];
T = 500;
metrics = zeros(length(factors) * length(sample_factors), T);
n_region = zeros(length(factors) * length(sample_factors), T);
parfor seed = 1:T
    [metrics(:, seed), n_region(:, seed)] = sim4_all_seeds(loc, radius, base_num_in_circle, factors, lambda, sample_factors, seed);
end
save('sim4_all_seeds_result.mat')
   
