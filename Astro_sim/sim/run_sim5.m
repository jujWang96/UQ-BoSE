addpath(genpath('/home/minjay/G-SRG'))
addpath(genpath('/home/minjay/Astro_sim'))

range_x = [0.2 0.4; 0.4 0.6; 0.6 0.8];
range_y = [0.2 0.4; 0.2 0.8; 0.6 0.8];
loc = [0.35 0.3; 0.5 0.5; 0.65 0.7];
radius = [0.025 0.025 0.025];

base_num_in_L_shape = [2 6 2];
base_num_in_circle = ones(1, 3);
factors = [10 20 30];
lambda = 1000;
sample_factors = [0.5 1 2];
T = 500;

metrics = zeros(length(factors) * length(sample_factors), T);
n_region = zeros(length(factors) * length(sample_factors), T);
parfor seed = 1:T
    [metrics(:, seed), n_region(:, seed)] = sim5(range_x, range_y, base_num_in_L_shape,...
        loc, radius, base_num_in_circle, factors, lambda, sample_factors, seed);
end
save('sim5_result.mat')