addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))
addpath(genpath('~/Desktop/gsrg/plot_util'))

loc = [0.5 0.5; 0.4 0.4; 0.4 0.6; 0.6 0.4; 0.6 0.6];
radius = [0.25 0.025*ones(1, 4)];
base_num_in_circle = [10 ones(1, 4)];
factors = [10 20 30];
lambda = 1000;
sample_factors = [0.5 1 2];
T = 500;
metrics_beam = zeros(length(factors) * length(sample_factors), T);
n_region_beam = zeros(length(factors) * length(sample_factors), T);
merge_time_beam = zeros(length(factors) * length(sample_factors), T);
seedsize = 5;
for gridsize = [0.2,0.1,0.05]
    for seed = 1:T
            [metrics_beam(:, seed), n_region_beam(:, seed), merge_time_beam(:,seed)] = sim4_beam(loc, radius, base_num_in_circle, factors, lambda, sample_factors, seed,gridsize, seedsize);
       
    end

save(strcat('sim4_result_beam_seedsize5_gs',num2str(gridsize),'.mat'))
end 
metrics_beam = zeros(length(factors) * length(sample_factors), T);
n_region_beam = zeros(length(factors) * length(sample_factors), T);
merge_time_beam = zeros(length(factors) * length(sample_factors), T);
for gridsize = [0.2,0.1,0.05]
    seedsize = ceil(25/(1/gridsize));
    parfor seed = 1:T
            [metrics_beam(:, seed), n_region_beam(:, seed), merge_time_beam(:,seed)] = sim4_beam(loc, radius, base_num_in_circle, factors, lambda, sample_factors, seed,gridsize, seedsize);
       
    end
end
save(strcat('sim4_result_beam_gs',num2str(gridsize),'.mat'))
   
