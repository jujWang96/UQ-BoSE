
% Circular extended source
clear
close all
addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))
addpath(genpath('~/Desktop/gsrg/plot_util'))
bw = 10;

GRAY = [0.6 0.6 0.6];

loc = [0.5 0.5; 0.4 0.4; 0.4 0.6; 0.6 0.4; 0.6 0.6];
radius = [0.25 0.025*ones(1, 4)];

base_num_in_circle = [10 ones(1, 4)];

factor = 30;
sample_factor = 1;
lambda = 10000;
seed = 9;

X = sim_inhomo_Pois_const([0 1], [0 1], lambda, loc, radius, factor * base_num_in_circle, seed);


subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);

%specify the line width of the circles
lw = 2;



% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[seeds, seeds_rej, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(0.1, 0.9, 0.1, 0.9,...
    0.2, 0.2, 5, cell_log_intensity, cell_area, cx, cy, 2, 50, 5, invalid, adj_mat);
num = num_s+num_s_pt;
disp(['Number of regions is ', num2str(num)])

step = num;
seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');
[region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);

tic
[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);
toc
BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);


tic
[sets_all, log_like_all] = merge_region_fast(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n);
toc
BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);

%specify the color map for background: gray, point: blue, extend: red



tic
[log_like_buffer,sets_buffer,sets_all, log_like_all,region_area_buffer,region_intensity_buffer, region_num_cells_buffer] = merge_region_beam_search(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n, bw,step);
toc
BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);

%specify the color map for background: gray, point: blue, extend: red



tic
[log_like_buffer_fb,sets_buffer_fb,sets_all_fb, log_like_all_fb,region_area_buffer_fb,region_intensity_buffer_fb, region_num_cells_buffer_fb] = merge_region_fast_beamsearch_optimization(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n, bw);
toc
BIC_all = -2*log_like_all_fb+4*(num-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);

%specify the color map for background: gray, point: blue, extend: red

debug_beamsearch