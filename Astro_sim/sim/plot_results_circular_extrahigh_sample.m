g1 = [ones(2000, 1); 2*ones(2000, 1); 3*ones(2000, 1); 4*ones(2000, 1)];
g2 = [ones(500, 4); 2*ones(500, 4); 3*ones(500, 4); 4*ones(500, 4)];
g2 = g2(:);
reduction = 0.02;
addpath(genpath('../sim_util'))
addpath(genpath('../util'))

% plot for sim4
load('sim4_result_extra_seedsize5_gs0.2.mat')
figure
%boxplot_group(metrics', g1, g2, '')
xlabel('\beta')
%min_white_margin(gca, -reduction/2, reduction)
%saveas(gca,'metrics_circular_new.pdf')

figure
true_n_source = 6; %5;
histogram_plot_group(n_region, true_n_source, factors, sample_factors)
calc_true_source_num(n_region, true_n_source, factors, sample_factors)
%saveas(gca,'nsource_circular_new.pdf')

