addpath(genpath('~/UQ-BoSE/Astro_sim'))
addpath(genpath('~/UQ-BoSE/Astro_sim/sim_util'))
addpath(genpath('~/UQ-BoSE/G-SRG'))
addpath(genpath('~/UQ-BoSE/plot_util'))

loc = [0.5 0.5; 0.4 0.4; 0.4 0.6; 0.6 0.4; 0.6 0.6];
radius = [0.25 0.025*ones(1, 4)];
base_num_in_circle = [10 ones(1, 4)];
factors = [10 20 30];
lambda = 1000;
sample_factors = [0.5 1 2];
T = 500;
metrics_random = zeros(length(factors) * length(sample_factors), T);
n_region_random = zeros(length(factors) * length(sample_factors), T);
merge_time_random = zeros(length(factors) * length(sample_factors), T);

for gridsize = [0.2, 0.1,0.05]
    seedsize = ceil(25/(1/gridsize));
    for seed = 1:T
       
        [metrics_random(:, seed), n_region_random(:, seed), merge_time_random(:,seed)] = sim4_random(loc, radius, base_num_in_circle, factors, lambda, sample_factors, seed,gridsize,seedsize);
        
    end
    save(strcat('sim4_result_random_gs',num2str(gridsize),'.mat'))
end

seedsize = 5;

for gridsize = [0.2,0.1,0.05]
    parfor seed = 1:T
       
        [metrics_random(:, seed), n_region_random(:, seed), merge_time_random(:,seed)] = sim4_random(loc, radius, base_num_in_circle, factors, lambda, sample_factors, seed,gridsize,seedsize);
        
    end
    save(strcat('sim4_result_random_seedsize5_gs',num2str(gridsize),'.mat'))
end       