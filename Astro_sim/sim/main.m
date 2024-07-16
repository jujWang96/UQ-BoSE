% Circular extended source

clear
close all

 addpath(genpath('../sim_util'))
 addpath(genpath('../util'))
 addpath(genpath('../../G-SRG-master'))



load notconnect.mat
X = X2;

% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[seeds, seeds_rej, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(0.1, 0.9, 0.1, 0.9,...
    0.05, 0.05, 5, cell_log_intensity, cell_area, cx, cy, 2, 20, 30, invalid,adj_mat);

num = num_s+num_s_pt;
disp(['Number of regions is ', num2str(num)])

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');
reg1 = region_sets{141};

GRAY = [0.6 0.6 0.6];
triplot(DT, 'Color', GRAY)
axis equal
hold on 
scatter(X(:,1),X(:,2),'.');
scatter(X(reg1,1),X(reg1,2),25,'k');
