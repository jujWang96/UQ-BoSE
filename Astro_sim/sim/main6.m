% Arc-shaped extended source
clear
close all
addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))
addpath(genpath('~/Desktop/gsrg/plot_util'))
addpath(genpath('~/Desktop/gsrg/contour_util'))
clear
close all

GRAY = [0.6 0.6 0.6];

loc_ring = [0.3 0.5];
radius_in = 0.2;
radius_out = 0.4;
loc = [0.5 0.7; 0.6 0.5; 0.5 0.3];
radius = 0.025*ones(1, 3);

base_num_in_ring = 10;
base_num_in_circle = ones(1, 3);

factor = 30;
sample_factor = 1;
lambda = 1000;
seed = 0;

X = sim_inhomo_Pois_const_ring(loc_ring, radius_out, radius_in, factor * base_num_in_ring, seed);
X = [X; sim_inhomo_Pois_const([0 1], [0 1], lambda, loc, radius, factor * base_num_in_circle)];
csvwrite('X_zshape.csv',X)
h = figure;

%specify the linewidth
lw_circle = 2;
lw_arc = 2;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);

subplot(1, 3, 1)
plot_circles(loc, radius,lw_circle)
hold on
plot_arc(loc_ring, radius_in, radius_out,lw_arc)
scatter(X(:, 1), X(:, 2), 'k.')
axis([0 1 0 1])
axis square
box on
set(gca, 'fontsize', 12)

% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[seeds, seeds_rej, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(0.1, 0.9, 0.1, 0.9,...
    0.2, 0.2, 5, cell_log_intensity, cell_area, cx, cy, 2, 50, 5, invalid, adj_mat);
num = num_s+num_s_pt;
disp(['Number of regions is ', num2str(num)])

% plot the seeds
subplot(1, 3, 2)
% specify the colormap
colors = lines(num);
plot_seeds2(DT, cx, cy, seeds, seeds_pt, seeds_rej, colors, num_s, num_s_pt)
set(gca, 'fontsize', 12)

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);

BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);

subplot(1, 3, 3)
%specify the color map for background: gray, point: blue, extend: red
backg_seg = 1;
extend_seg = 7;
point_seg = [12,14,39];
colors = [GRAY;[0.8500,0.3250,0.0980];repmat([0,0.4470,0.7410],3,1)];

plot_segmentation(DT, index_BIC, sets_all, cx, cy, colors)
hold on
plot_circles(loc, radius,lw_circle)
plot_arc(loc_ring, radius_in, radius_out,lw_arc)
axis image
set(gca, 'fontsize', 12)

set(h, 'Position', [0, 0, 800, 260]);
