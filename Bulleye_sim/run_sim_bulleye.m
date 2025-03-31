addpath(genpath('~/UQ-BoSE/Bulleye_sim'))
addpath(genpath('~/UQ-BoSE/Astro_sim'))
addpath(genpath('~/UQ-BoSE/G-SRG'))
addpath(genpath('~/UQ-BoSE/plot_util'))
addpath(genpath('~/UQ-BoSE/contour_util'))
addpath(genpath('~/UQ-BoSE/FD_util'))


dist = 0;
loc = [0.5+sign(dist)*sqrt(dist^2/2) 0.5+sign(dist)*sqrt(dist^2/2)];
unit_a = 0.03;
%radius = [(unit_a/pi)^(1/2),(4*unit_a/pi)^(1/2)-(unit_a/pi)^(1/2),(9*unit_a/pi)^(1/2)-(4*unit_a/pi)^(1/2)];
radius = [(unit_a/pi)^(1/2),(4*unit_a/pi)^(1/2)-(unit_a/pi)^(1/2)];

lambdas = [100, 500,2500];
 contrasts = [1.4,2,2.8];
gridszie = 0.15;
seedsize = 1;
T = 500;
rep_itr = 1000;
SNR = 3;
metrics = zeros(length(lambdas) * length(contrasts), T);
n_region = zeros(length(lambdas) * length(contrasts), T);
rad = zeros(length(lambdas) * length(contrasts), T);
parfor seed = 1:T
    [metrics(:,seed), n_region(:,seed),rad(:,seed)]= sim_bulleye(loc, radius, unit_a, SNR, lambdas, contrasts, seed,gridszie,seedsize,rep_itr);  
end

save(strcat('sim_bulleye_result_SNR',num2str(SNR),'.mat'))

SNR = 10;
metrics = zeros(length(lambdas) * length(contrasts), T);
n_region = zeros(length(lambdas) * length(contrasts), T);
rad = zeros(length(lambdas) * length(contrasts), T);
parfor seed = 1:T
    [metrics(:,seed), n_region(:,seed),rad(:,seed)]= sim_bulleye(loc, radius, unit_a, SNR, lambdas, contrasts, seed,gridszie,seedsize,rep_itr);  
end

save(strcat('sim_bulleye_result_SNR',num2str(SNR),'.mat'))

SNR = 30;
metrics = zeros(length(lambdas) * length(contrasts), T);
n_region = zeros(length(lambdas) * length(contrasts), T);
rad = zeros(length(lambdas) * length(contrasts), T);
parfor seed = 1:T
    [metrics(:,seed), n_region(:,seed),rad(:,seed)]= sim_bulleye(loc, radius, unit_a, SNR, lambdas, contrasts, seed,gridszie,seedsize,rep_itr);  
end

save(strcat('sim_bulleye_result_SNR',num2str(SNR),'.mat'))
