clear
close all
load('real_2_voronoi_result_penalty6_2023_2_12_22_16_2.mat')
addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))
addpath(genpath('~/Desktop/gsrg/plot_util'))
filename = 'real_2_random_voronoi';
filename = strcat(filename,'_penalty',num2str(penalty_term));

% plot the data
fig = figure;
scatter(X(:, 1), X(:, 2), 1,'k.')
axis image
box on
set(gca, 'fontsize', 18)
min_white_margin(gca);
set(fig, 'Position', [0, 0, 490, 380]);

saveas(fig, strcat(filename,'data'), 'epsc')

adj_mat = get_adj_mat( E, n );
% plot the log intensity
fig = figure;
triplot(DT, 'Color', GRAY)
hold on

[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
on_margin = find(cell_area> mean(cell_area(valid))+5*(var(cell_area(valid)).^(1/2)));
valid = setdiff(valid,on_margin);
invalid = [invalid; on_margin];

validpt = [cx(valid),cy(valid),cell_log_intensity(valid)];
validpt = sortrows(validpt,3);
% scatter(cx(valid), cy(valid), 20, cell_log_intensity(valid), 'filled')
scatter(validpt(:,1), validpt(:,2), 20, validpt(:,3), 'filled')

%c.Location('east')
colormap(jet)
c = colorbar;
c.Location = 'manual';
axis image
set(gca, 'fontsize', 18)
min_white_margin(gca);
set(fig, 'Position', [0, 0, 490, 380]);
set(c,'Position',[0.924,  0.073, 0.03, 0.9])


saveas(fig, strcat(filename,'log_intensity_in'), 'epsc')

% plot the seeds
fig = figure;
colors = lines(num);
plot_seeds_one(DT, cx, cy, seeds, [], [], colors, num, 0)
set(gca, 'fontsize', 18)
min_white_margin(gca);
set(fig, 'Position', [0, 0, 490, 380]);

saveas(fig, strcat(filename,'seeds'), 'epsc')

% drop empty entries in the cell array
num = length(selected);

area_all = [];
log_int_all = [];
selected_nonempty = {};
index = 0;
for i = 1:num
    if ~isempty(selected{i})
        index = index+1;
        selected_nonempty{index} = selected{i};
        area = sum(cell_area(selected{i}));
        area_all = [area_all, area];
        log_int = log(sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i}))/area);
        log_int_all = [log_int_all, log_int];
    end
end

num_nonempty = length(selected_nonempty);

% get source boundaries formed by Voronoi cell edges
[V, R] = voronoiDiagram(DT);
vx_edges_all = {};
vy_edges_all = {};
for j = 1:num_nonempty
    cells_in_region = R(selected_nonempty{j});
    edges = [];
    for i = 1:length(cells_in_region)
        % each row records two vertex indices of an edge
        new_edges = [cells_in_region{i}' [cells_in_region{i}(2:end)'; cells_in_region{i}(1)]];
        % sort vertex indices to avoid ambiguity
        new_edges = sort(new_edges, 2);
        edges = [edges; new_edges];
    end
    [unique_edges, ~ , ind] = unique(edges, 'rows');
    % get the count of each unique edge
    counts = histc(ind, unique(ind));
    % remove edges that appear at least twice to get the boundary edges
    edges = unique_edges(counts==1, :);
    % x-axis
    vx_edges = [V(edges(:, 1), 1)'; V(edges(:, 2), 1)'];
    % y-axis
    vy_edges = [V(edges(:, 1), 2)'; V(edges(:, 2), 2)'];
    nume = size(vx_edges, 2);
    vx_edges = [vx_edges; NaN(1, nume)];
    vx_edges = vx_edges(:);
    vx_edges_all{j} = vx_edges;
    vy_edges = [vy_edges; NaN(1, nume)];
    vy_edges = vy_edges(:);
    vy_edges_all{j} = vy_edges;
end

fig = figure;
scatter(cx, cy, 3, repmat(GRAY, length(cx), 1), 'filled')

hold on
scatter(X(invalid,1),X(invalid,2))
for j = 1:num_nonempty
    % original red is brighter
    line(vx_edges_all{j}, vy_edges_all{j}, 'Color', 'r')
end
axis image
set(gca, 'fontsize', 18)
box on
min_white_margin(gca);
saveas(fig, strcat(filename,'segment_result'), 'epsc')


left_corner1 = [0.2 0.45];
right_corner1 = [0.45 0.6];
rectangle('Position', [left_corner1 right_corner1 - left_corner1], 'LineWidth', 1)
text(left_corner1(1), right_corner1(2)+0.02, 'Region 1', 'FontSize', 16)

left_corner2 = [0.6 0.3];
right_corner2 = [0.8 0.5];
rectangle('Position', [left_corner2 right_corner2 - left_corner2], 'LineWidth', 1)
text(left_corner2(1), right_corner2(2)+0.02, 'Region 2', 'FontSize', 16)

left_corner3 = [0.3 0.05];
right_corner3 = [0.55 0.2];
rectangle('Position', [left_corner3 right_corner3 - left_corner3], 'LineWidth', 1)
text(left_corner3(1), right_corner3(2)+0.02, 'Region 3', 'FontSize', 16)


fig = figure;
scatter(cx, cy, 10, repmat(GRAY, length(cx), 1), 'filled')
hold on
for j = 1:num_nonempty
    line(vx_edges_all{j}, vy_edges_all{j}, 'Color', 'r', 'LineWidth', 1)
end
axis image
axis([left_corner1(1) right_corner1(1) left_corner1(2) right_corner1(2)])
title('Region 1')
set(gca, 'FontSize', 18)
box on
%min_white_margin(gca);

set(fig, 'Position', [0, 0, 425, 400]);

saveas(fig, strcat(filename,'region1.png') )

fig = figure;
scatter(cx, cy, 10, repmat(GRAY, length(cx), 1), 'filled')
hold on
for j = 1:num_nonempty
    line(vx_edges_all{j}, vy_edges_all{j}, 'Color', 'r', 'LineWidth', 1)
end
axis image
axis([left_corner2(1) right_corner2(1) left_corner2(2) right_corner2(2)])
title('Region 2')
set(gca, 'FontSize', 18)
box on
%min_white_margin(gca);
set(fig, 'Position', [0, 0, 425, 400]);
saveas(fig, strcat(filename,'region2.png'))

fig = figure;
scatter(cx, cy, 10, repmat(GRAY, length(cx), 1), 'filled')
hold on
for j = 1:num_nonempty
    line(vx_edges_all{j}, vy_edges_all{j}, 'Color', 'r', 'LineWidth', 1)
end
axis image
axis([left_corner3(1) right_corner3(1) left_corner3(2) right_corner3(2)])
title('Region 3')
set(gca, 'FontSize', 18)
box on
%min_white_margin(gca);
set(fig, 'Position', [0, 0, 425, 400]);
saveas(fig, strcat(filename,'region3.png'))


fig = figure;
subplot(3,1,1)
scatter(cx, cy, 10, repmat(GRAY, length(cx), 1), 'filled')
hold on
for j = 1:num_nonempty
    line(vx_edges_all{j}, vy_edges_all{j}, 'Color', 'r', 'LineWidth', 1)
end
axis image
axis([left_corner1(1) right_corner1(1) left_corner1(2) right_corner1(2)])
title('Region 1')
set(gca, 'FontSize', 18)
box on
min_white_margin(gca);
subplot(3,1,2)

scatter(cx, cy, 10, repmat(GRAY, length(cx), 1), 'filled')
hold on
for j = 1:num_nonempty
    line(vx_edges_all{j}, vy_edges_all{j}, 'Color', 'r', 'LineWidth', 1)
end
axis image
axis([left_corner2(1) right_corner2(1) left_corner2(2) right_corner2(2)])
title('Region 2')
set(gca, 'FontSize', 18)
box on
set(fig, 'Position', [0, 0, 425, 400]);
min_white_margin(gca);
subplot(3,1,3)

scatter(cx, cy, 10, repmat(GRAY, length(cx), 1), 'filled')
hold on
for j = 1:num_nonempty
    line(vx_edges_all{j}, vy_edges_all{j}, 'Color', 'r', 'LineWidth', 1)
end
axis image
axis([left_corner3(1) right_corner3(1) left_corner3(2) right_corner3(2)])
title('Region 3')
set(gca, 'FontSize', 18)
box on
min_white_margin(gca);
saveas(fig, strcat(filename,'region_stack'), 'epsc')

% point sources (to be compared with the result from the wavdetect
% algorithm)
% not to be fancy, choose a cutoff point
cutoff_point_source = 0.0003;

fig = figure;
triplot(DT, 'Color', GRAY)
hold on
for i = 1:num_nonempty
    % if area is less than the threshold
    if area_all(i) < cutoff_point_source
        scatter(cx(selected_nonempty{i}), cy(selected_nonempty{i}), 12, log_int_all(i)*ones(length(selected_nonempty{i}), 1), 'filled')
    end
end
colorbar('EastOutside')
colormap(hsv)
axis image
set(gca, 'fontsize', 14)
min_white_margin(gca);
saveas(fig, strcat(filename,'point_sources'), 'epsc')

% fig = figure;
% plot(num:-1:1, BIC_all + log(n), '-o', 'MarkerSize', 3)
% xlabel('Number of segments')
% ylabel('BIC')
% y_range = get(gca, 'ylim');
% hold on
% plot([num-index_BIC+1 num-index_BIC+1],y_range, 'LineWidth', 1.5)
% axis tight
% position = get(gcf, 'Position');
% % shrink width and height
% set(gcf, 'Position', [position(1), position(2), position(3)/1.5, position(4)/1.5])
% min_white_margin(gca);
fig = figure;
plot_segmentation_wo_voronoi(DT, selected, cx, cy,colors ,8,true)
set(gca, 'fontsize', 14);
min_white_margin(gca);
saveas(fig, strcat(filename,'segmentat_result_region'), 'epsc')

fig = figure;
triplot(DT, 'Color', GRAY)
hold on
% the final result
%selected = selected_rep{index_rep};

index = 0;
selected_nonempty = {};
for i = 1:length(selected)
    if ~isempty(selected{i})
        index = index + 1;
        selected_nonempty{index} = selected{i};
    end
end

for i = 1:length(selected_nonempty)-1
    for j = i+1:length(selected_nonempty)
        % sorted by area from the largest to the smallest
        if sum(cell_area(selected_nonempty{i})) < sum(cell_area(selected_nonempty{j}))
            tmp = selected_nonempty{i};
            selected_nonempty{i} = selected_nonempty{j};
            selected_nonempty{j} = tmp;
        end
    end
end

for i = 1:length(selected_nonempty)
    log_int = log(sum(exp(cell_log_intensity(selected_nonempty{i})).*cell_area(selected_nonempty{i}))/sum(cell_area(selected_nonempty{i})));
    scatter(cx(selected_nonempty{i}), cy(selected_nonempty{i}), 12, log_int*ones(length(selected_nonempty{i}), 1), 'filled')
end
colorbar('EastOutside')
colormap(hsv)
caxis([9 16])
axis image
set(gca, 'fontsize', 12)
min_white_margin(gca);
saveas(fig, strcat(filename,'segment_result_points'), 'epsc')
