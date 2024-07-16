GRAY = [0.6 0.6 0.6];

[cy, cx] = find(y>0);
X = [cx cy];
bound_x = [min(cx) max(cx)];
bound_y = [min(cy) max(cy)];
count = y(y>0);

% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, bound_x, bound_y, count);

figure
triplot(DT, 'Color', GRAY)
hold on
invalid = find(isnan(cell_log_intensity));
valid = setdiff(1:n, invalid);
scatter(cx(valid), cy(valid), 12, cell_log_intensity(valid), 'filled')
colorbar
colormap(jet)
axis image

% get seeds
[seeds, seeds_rej, seeds_pt, num_s, num_s_pt, invalid] = get_seeds_sim_local_max(min(cx)+5, max(cx)-5, min(cy)+5, max(cy)-5,...
    100, 100, 5, cell_log_intensity, cell_area, cx, cy, 2, 100, 5);
num = num_s+num_s_pt;
disp(['Number of regions is ', num2str(num)])

% plot the seeds
figure
triplot(DT, 'Color', GRAY)
hold on
% specify the colormap
colors = lines(num_s);
for i = 1:num_s
    scatter(cx(seeds{i}), cy(seeds{i}), 12, colors(i, :), 's', 'filled')
end
for i = 1:num_s_pt
    scatter(cx(seeds_pt{i}), cy(seeds_pt{i}), 12, 'r', 'd', 'filled')
end
for i = 1:length(seeds_rej)
    scatter(cx(seeds_rej{i}), cy(seeds_rej{i}), 12, 'k', '^', 'filled')
end
axis image

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
adj_mat = get_adj_mat( E, n );
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

% plot the over-segmented image
figure
triplot(DT, 'Color', GRAY)
hold on
for i = 1:num_s
    scatter(cx(region_sets{i}), cy(region_sets{i}), 12,  colors(i, :), 'filled')
end
for i = 1:num_s_pt
    scatter(cx(region_sets{i+num_s}), cy(region_sets{i+num_s}), 12, 'r', 'filled')
end
axis image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);

%figure
%plot(log_like_all, '-o')

BIC_all = -2*log_like_all+6*(num-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);

figure
plot(num-1:-1:0, BIC_all, '-o', 'MarkerSize', 3)

figure
triplot(DT, 'Color', GRAY)
hold on
% the final result
selected = sets_all{index_BIC};
for i = 1:num
    if ~isempty(selected{i})
        log_int = log(sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i}))/sum(cell_area(selected{i})));
        scatter(cx(selected{i}), cy(selected{i}), 12, log_int*ones(length(selected{i}), 1), 'filled')
    end
end
colorbar('SouthOutside')
colormap(hsv)
axis image
