

function [labeled_cells, pred_class_all, n_region,merge_time] = sim_fit_gridonly(X, factor, sample_factor, seed, show_plot,gridszie,seedsize)
% fit the model in simulation

% init comp
if nargin == 6
    seedsize = ceil(25/(1/gridsize));
end
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

% get seeds
[invalid, ~] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[seeds, ~,  num_s] = get_seeds_sim(0.1, 0.9, 0.1, 0.9,...
       gridsize,gridsize, seedsize, cell_log_intensity, cell_area, cx, cy, 2, invalid,adj_mat);
num = num_s;
disp(['Number of regions is ', num2str(num)])


% make a copy of variable seeds
region_sets = seeds;


% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

if length(labeled_cells) + length(invalid) < n
    dropped_cells = setdiff(setdiff(1:n, labeled_cells), invalid);
    disp(['Premature termination happened at factor=', num2str(factor), ',sample_factor=', num2str(sample_factor), ',seed=', num2str(seed)])
    fig = figure;
    plot(cx, cy, '.')
    hold on
    plot(cx(dropped_cells), cy(dropped_cells), 'ro')
    saveas(fig, ['debug_', num2str(factor), '_', num2str(sample_factor), '_', num2str(seed)], 'png')
    save(['debug_', num2str(factor), '_', num2str(sample_factor), '_', num2str(seed), '.mat'])
end
tic;
[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);

% the factor is 4 instead of 6
BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
[~, index_BIC] = min(BIC_all);
merge_time = toc;

selected = sets_all{index_BIC};
n_region = 0;
pred_class_all = zeros(n, 1);
for i = 1:num
    if ~isempty(selected{i})
        n_region = n_region+1;
        pred_class_all(selected{i}) = n_region;
    end
end

if show_plot
    figure
    GRAY = [0.6 0.6 0.6];
    colors = lines(n_region);
    triplot(DT, 'Color', GRAY)
    hold on
    index = 0;
    for i = 1:num
        if ~isempty(selected{i})
            index = index+1;
            scatter(cx(selected{i}), cy(selected{i}), 12, colors(index, :), 'filled')
        end
    end
    axis image
end

end
