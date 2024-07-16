
close all
clear
addpath(genpath('~/src/gsrg/NGC2300'))
load('ngc2300_box_058kev_evt2.mat')

%X = double(unique(X,'rows'));
n = length(X);
% fraction = 1
% id = randsample(n,round(n*fraction));
% X = X(id,:);
X = double(unique(X,'rows'));
GRAY = [0.6,0.6,0.6];
P= 14;
threshold = 15;
M = 150;
imagename = '~/src/gsrg/Astro_sim/real_img/sim_results/NGC2300chandra_vor_vk_';
%X = double(unique(X,'rows'));
n = length(X);
% fraction = 1
% id = randsample(n,round(n*fraction));
% X = X(id,:);
X = double(unique(X,'rows'));
GRAY = [0.6,0.6,0.6];
%P= 17;
%threshold = 5;
%M = 150;

penalty = 6;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[V, R] = voronoiDiagram(DT);

% get seeds
[seeds, num_region] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat,cell_area,P, threshold);

seeds_all = seeds ;

% make a copy of variable seeds

disp(strcat('The number of seeds is',num2str(num_region)))


colors = lines(num_region);
figure
triplot(DT, 'Color', GRAY)
hold
for i = 1:num_region
    scatter(cx(seeds{i}), cy(seeds{i}), 12, colors(i, :), 's', 'filled')
end
axis equal

region_sets = seeds;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

[region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num_region, cell_area, cell_log_intensity, region_sets, adj_mat);
sets_greedy = region_sets;

%greedy merge until M region left then do randomized search
[sets_all_greedy, log_like_all_greedy]  = merge_region_fast(num_region, region_area, ...
region_intensity, region_sets, adj_mat_region, region_num_cells, n);

fig = figure;
n_region = num_region:-1:2;
d_log_like_all = diff(log_like_all_greedy);
plot(n_region, (d_log_like_all), '-o', 'MarkerSize', 3)
xlabel('Number of segments')
ylabel('\delta loglikelihood')

hold on
axis tight
position = get(gcf, 'Position');
% shrink width and height
set(gcf, 'Position', [position(1), position(2), position(3)/1.5, position(4)/1.5])
min_white_margin(gca);
saveas(gcf, strcat(imagename, 'delta_loglikelihood'),'epsc')

fig = figure;
n_region = num_region:-1:2;
d_log_like_all = diff(log_like_all_greedy);
plot(n_region, exp(d_log_like_all), '-o', 'MarkerSize', 3)
xlabel('Number of segments')
ylabel('exp(\delta loglikelihood)')

hold on
axis tight
position = get(gcf, 'Position');
% shrink width and height
set(gcf, 'Position', [position(1), position(2), position(3)/1.5, position(4)/1.5])
min_white_margin(gca);
saveas(gcf, strcat(imagename, 'delta_loglikelihood'),'epsc')
