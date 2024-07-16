% numerical experiment 1
% one round-shaped extended source

clear
close all

% generate simulated data (inhomogeneous Poisson point process)
X = sim_inhomo_Pois_const([0 1], [0 1], 1000, [0.5 0.5], 0.25, 200);

% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1]);

% plot log intensity
figure
triplot(DT)
hold on
scatter(cx, cy, [], cell_log_intensity, 'o', 'filled')
colorbar
colormap(jet)
axis image

% get seeds
[seeds, num, invalid] = get_seeds_sim(0.1, 0.9, 0.1, 0.9,...
    0.2, 0.2, 5, cell_log_intensity, cell_area, cx, cy, 2);
disp(['Number of regions is ', num2str(num)])

% plot the seeds
figure
triplot(DT)
hold on
% specify the colormap
colors = lines(num);
for i = 1:num
    scatter(cx(seeds{i}), cy(seeds{i}), [], colors(i, :), '*')
end
axis image

% make a copy of variable seeds
region_sets = seeds;

% graph-based SRG
adj_mat = get_adj_mat( E, n );
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

% plot the over-segmented image
figure
triplot(DT)
hold on
for i = 1:num
    scatter(cx(region_sets{i}), cy(region_sets{i}), [],  colors(i, :), 'filled')
end
viscircles([0.5 0.5], 0.25, 'EdgeColor', 'r', 'LineWidth', 1.5);
axis image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);

figure
plot(log_like_all, '-o')

BIC_all = -2*log_like_all+6*(num-1:-1:0)'*log(n);
figure
plot(num-1:-1:0, BIC_all, '-o')
xlabel('Number of clusters/sources')
ylabel('BIC')

[~, index_BIC] = min(BIC_all);

figure
triplot(DT)
hold on
% the final result
selected = sets_all{index_BIC};
index = 0;
for i = 1:num
    if ~isempty(selected{i})
        index = index+1;
        scatter(cx(selected{i}), cy(selected{i}), [],  colors(index, :), 'filled')
    end
end
viscircles([0.5 0.5], 0.25, 'EdgeColor', 'r', 'LineWidth', 1.5);
axis image