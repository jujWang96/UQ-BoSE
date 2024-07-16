

function [labeled_cells, pred_class_all, n_region,merge_time] = sim_fit_random_gridonly(X, factor, sample_factor, seed, show_plot,gridsize,seedsize)
% fit the model in simulation with random search 
rng(seed)
rand_num = 3;
rep_itr = 100;
if nargin == 6
    seedsize = ceil(25/(1/gridsize));
end
% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

% get seeds
[invalid, ~] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[seeds, ~,  num_s] = get_seeds_sim(0.1, 0.9, 0.1, 0.9,...
       gridsize,gridsize, seedsize, cell_log_intensity, cell_area, cx, cy, 2, invalid,adj_mat);
 num = num_s;
%disp(['Number of regions is ', num2str(num)])

seeds_all = seeds ;

% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

[region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);

selected_rep = cell(1,rep_itr);
BIC_rep = zeros(1,rep_itr);
tic;
%repeat multiple times
parfor ii = 1:rep_itr
    [sets_all, log_like_all] = merge_region_random_fast(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n, rand_num);
    % the factor is 4 instead of 6
    BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
    [min_BIC, index_BIC] = min(BIC_all);
    selected = sets_all{index_BIC};
    selected_rep{ii} = selected;
    BIC_rep(ii) = min_BIC;
end

[min_BIC, index_rep] = min(BIC_rep);
selected = selected_rep{index_rep};
merge_time = toc;
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
