function [impute_log_int, n_region] = real_full_jitter(seed)
rng(seed)

load('real_full_result_2018_11_1_23_4_37.mat')

cx_origin = cx;
cy_origin = cy;

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

% jitter X
[sx_all, sy_all, ~] = get_bootstrap_samples(cx, cy,...
    DT, num_nonempty, selected_nonempty, log_int_all);
X_jitter = [sx_all', sy_all'];
% drop duplicate rows
X_jitter = unique(X_jitter, 'rows');
n = size(X_jitter, 1);
count = ones(n, 1);
[cx, cy, n, ~, E, cell_log_intensity, cell_area] = init_comp(X_jitter, bound_x, bound_y, count);
adj_mat = get_adj_mat( E, n );

[invalid, ~] = get_invalid_cells(cell_log_intensity, adj_mat, n);

disp('Getting initial seeds...')
% get seeds
st_x = 0.1;
en_x = 0.9;
st_y = 0.1*bound_y(2);
en_y = 0.9*bound_y(2);
step_x = 0.1;
step_y = 0.1*bound_y(2);
set_size = 20;
factor = 2;
k = 100;
set_size2 = 20;
[seeds, ~, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(st_x, en_x, st_y, en_y,...
    step_x, step_y, set_size, cell_log_intensity, cell_area, cx, cy, factor, k, set_size2, invalid);
num = num_s+num_s_pt;
disp(['Number of regions is ', num2str(num)])

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

disp('Seeded region growing...')
% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid', true, 1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Region merging...')
[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);

BIC_all = -2*log_like_all+6*(num-1:-1:0)'*log(n);
[~, index_BIC] = min(BIC_all);

cx_labeled = cx(labeled_cells);
cy_labeled = cy(labeled_cells);

selected = sets_all{index_BIC};
num = length(selected);
% Log intensity of all resampled Voronoi cells.
log_int_all = zeros(length(cx), 1);
n_region = 0;
for i = 1:num
    if ~isempty(selected{i})
        n_region = n_region + 1;
        area = sum(cell_area(selected{i}));
        log_int = log(sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i}))/area);
        log_int_all(selected{i}) = log_int;
    end
end

% Impute log intensity of all original Voronoi cells. 
impute_log_int = zeros(length(cx_origin), 1);
for i = 1:length(cx_origin)
    % Find the nearest (in terms of Euclidean distance) valid resampled
    % Voronoi cell to impute each original Voronoi cell.
    dist_vec = pdist2([cx_origin(i), cy_origin(i)], [cx_labeled, cy_labeled]);
    [~, index] = min(dist_vec);
    impute_log_int(i) = log_int_all(labeled_cells(index));
end

end
