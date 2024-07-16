addpath(genpath('/home/minjay/G-SRG'))
addpath(genpath('/home/minjay/Astro_sim'))

loc_ring = [0.3 0.5];
radius_out = 0.4;
radius_in = 0.2;
base_num_in_ring = 10;

loc = [0.5 0.7; 0.6 0.5; 0.5 0.3];
radius = 0.025*ones(1, 3);
base_num_in_circle = ones(1, 3);
lambda = 1000;
factors = [10 20 30];
sample_factors = [0.5 1 2];
T = 500;

metrics = zeros(length(factors) * length(sample_factors), T);
n_region = zeros(length(factors) * length(sample_factors), T);
parfor seed = 1:T
    [metrics(:, seed), n_region(:, seed)] = sim6(loc_ring, radius_out, radius_in,...
        base_num_in_ring, loc, radius, base_num_in_circle, lambda, factors,...
        sample_factors, seed);
end
save('sim6_result.mat')