tic
addpath(genpath('~/Desktop/gsrg/G-SRG'))
addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/plot_util'))

GRAY = [0.6 0.6 0.6];

% set global random seed
rng(0)

% read data
filename='photon_loc.txt';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

X = A.data;
clear A

% normalize locations
min_x = min(X(:, 1));
max_x = max(X(:, 1));
min_y = min(X(:, 2));
max_y = max(X(:, 2));
% the denominators have to be the same; otherwise, the area is distorted
X(:, 1) = (X(:, 1) - min_x) ./ (max_x - min_x);
X(:, 2) = (X(:, 2) - min_y) ./ (max_x - min_x);

disp('Conducting some initial computations...')
% init comp
bound_x = [0 1];
bound_y = [0 max(X(:, 2))];
n = size(X, 1);
count = ones(n, 1);
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, bound_x, bound_y, count);
adj_mat = get_adj_mat( E, n );

[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
region_sets = num2cell(valid);
region_area = cell_area(valid);
region_intensity = exp(cell_log_intensity(valid));
num_region = length(region_sets);
region_num_cells = ones(num_region, 1);

disp('Region merging...')
% TODO: Only return the best segmentation result instead of all
% segmentation history.
[sets_all, log_like_all] = merge_region_fast(num_region, region_area, ...
    region_intensity, region_sets, adj_mat(valid, valid), region_num_cells, n);

BIC_all = -2*log_like_all+6*(num_region-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);

fig = figure;
plot(num_region-1:-1:0, BIC_all, '-o', 'MarkerSize', 3)
saveas(fig, 'BIC', 'png')

fig = figure;
triplot(DT, 'Color', GRAY)
hold on
% the final result
selected = sets_all{index_BIC};
for i = 1:num_region
    if ~isempty(selected{i})
        log_int = log(sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i}))/sum(cell_area(selected{i})));
        scatter(cx(selected{i}), cy(selected{i}), 12, log_int*ones(length(selected{i}), 1), 'filled')
    end
end
colorbar('SouthOutside')
colormap(hsv)
axis image
saveas(fig, 'final_seg', 'png')

elapsed_time = toc;
disp(['Elapsed time is ', num2str(elapsed_time), ' seconds.'])

% save all variables (exclude figure handle)
clear fig
% it is too long to save
clear sets_all
filename = 'real_full_all_seeds_result';
current_datetime = clock;
for i = current_datetime
    filename = [filename, '_', num2str(fix(i))];
end
filename = [filename, '.mat'];
save(filename)
