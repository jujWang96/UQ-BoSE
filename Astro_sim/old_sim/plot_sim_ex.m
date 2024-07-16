% plot simulation examples
clear
close all

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.01], [0.05 0.05], [0.03 0.01]);

% record current number of plots
num_of_plot = 0;
titles = {'Low', 'Medium', 'High'};

% experiment 1
lambda_num = [1000 200; 1000 400; 1000 600];
loc = [0.5 0.5];
radius = 0.25;
num_of_plot = num_of_plot+1;
subplot(3, 4, num_of_plot)
viscircles(loc, radius, 'EdgeColor', 'k', 'LineWidth', 1);
axis([0 1 0 1])
axis square
title('True')
for i = 1:size(lambda_num, 1)
    lambda = lambda_num(i, 1);
    num_of_photon = lambda_num(i, 2);
    num_of_plot = num_of_plot+1;
    subplot(3, 4, num_of_plot)
    X = sim_inhomo_Pois_const([0 1], [0 1], lambda, loc, radius, num_of_photon, 1);
    scatter(X(:, 1), X(:, 2), '.')
    axis([0 1 0 1])
    axis square
    title(titles{i});
end

% experiment 2
lambda_num = [2000 50; 2000 100; 2000 150];
loc = [];
for i = 1:4
    for j = 1:4
        loc = [loc; i*0.2 j*0.2];
    end
end
radius = 0.05*ones(1, 16);
num_of_plot = num_of_plot+1;
subplot(3, 4, num_of_plot)
for i = 1:length(radius)
    viscircles(loc(i, :), radius(i), 'EdgeColor', 'k', 'LineWidth', 1);
end
axis([0 1 0 1])
axis square
for i = 1:size(lambda_num, 1)
    lambda = lambda_num(i, 1);
    num_of_photon_one = lambda_num(i, 2);
    num_of_photon = num_of_photon_one*ones(1, 16);
    num_of_plot = num_of_plot+1;
    subplot(3, 4, num_of_plot)
    X = sim_inhomo_Pois_const([0 1], [0 1], lambda, loc, radius, num_of_photon, 1);
    scatter(X(:, 1), X(:, 2), '.')
    axis([0 1 0 1])
    axis square
end

% experiment 3
lambda_num = [1000 100; 1000 200; 1000 300];
loc = [0.3 0.3; 0.5 0.7; 0.7 0.3];
radius = [0.15 0.1 0.05];
num_of_plot = num_of_plot+1;
subplot(3, 4, num_of_plot)
for i = 1:length(radius)
    viscircles(loc(i, :), radius(i), 'EdgeColor', 'k', 'LineWidth', 1);
end
axis([0 1 0 1])
axis square
for i = 1:size(lambda_num, 1)
    lambda = lambda_num(i, 1);
    num_of_photon_one = lambda_num(i, 2);
    num_of_photon = num_of_photon_one*ones(1, 3);
    num_of_plot = num_of_plot+1;
    subplot(3, 4, num_of_plot)
    X = sim_inhomo_Pois_const([0 1], [0 1], lambda, loc, radius, num_of_photon, 1);
    scatter(X(:, 1), X(:, 2), '.')
    axis([0 1 0 1])
    axis square
end
