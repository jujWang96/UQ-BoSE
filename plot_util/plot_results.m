g1 = [ones(1500, 1); 2*ones(1500, 1); 3*ones(1500, 1)];
g2 = [ones(500, 3); 2*ones(500, 3); 3*ones(500, 3)];
g2 = g2(:);
reduction = 0.02;

% plot for sim4
load('sim4_result.mat')
figure
boxplot_group(metrics', g1, g2, '')
xlabel('\alpha')
min_white_margin(gca, -reduction/2, reduction)

figure
true_n_source = 5;
histogram_plot_group(n_region, true_n_source, factors, sample_factors)

% plot for sim5
load('sim5_result.mat')
figure
boxplot_group(metrics', g1, g2, '')
xlabel('\alpha')
min_white_margin(gca, -reduction/2, reduction)

figure
true_n_source = 4;
histogram_plot_group(n_region, true_n_source, factors, sample_factors)

% plot for sim6
load('sim6_result.mat')
figure
boxplot_group(metrics', g1, g2, '')
xlabel('\alpha')
min_white_margin(gca, -reduction/2, reduction)

figure
true_n_source = 4;
histogram_plot_group(n_region, true_n_source, factors, sample_factors)
