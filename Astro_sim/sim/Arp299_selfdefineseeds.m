formatSpec='%f %f %d %d %d\n';
fileID = fopen('seeds_arp99.txt','r');
A = fscanf(fileID,formatSpec);
A = reshape(A,5,[])';
load('Arp299_merged_evt2_box_unscaled.mat')
load('Arp299_merged_evt2_box.mat')
addpath(genpath('~/Desktop/experiment/G-SRG-master'))
addpath(genpath('../sim'))
addpath(genpath('../sim_util'))
addpath(genpath('../util'))

X = double(unique(X,'rows'));

X_unscaled = double(unique(X_unscaled,'rows'));

[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);


%%%get the photon id of the seeds for all 
seeds_all = {};
for i = 1:length(A)
    [tf, index]=ismember(round(X_unscaled,2),A(i,1:2),'rows');
    %do not initialize the seeds if the photon is invalid
    if ~isempty(find(tf, 1))&&any(invalid==find(tf, 1))
        continue;
    end
    seeds_all{i} = find(tf, 1);
end
seeds_all = seeds_all(~cellfun('isempty',seeds_all));

num = length(seeds_all);
% get seeds

% plot the seeds
subplot(1, 2, 1)
% specify the colormap
colors = lines(num);
plot_seeds(DT, cx, cy, seeds_all, [], [], colors, num, 0)
set(gca, 'fontsize', 12)


% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);

BIC_all = -2*log_like_all+10*(num-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);

subplot(1, 3, 3)
plot_segmentation(DT, selected, cx, cy, lines(length(cx)))
hold on
plot_circles(loc, radius)
axis image
set(gca, 'fontsize', 12)

set(h, 'Position', [0, 0, 800, 260]);
