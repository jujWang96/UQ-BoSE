function [labeled_cells, pred_class_all, n_region] = sim_fit_all_seeds(X)
% fit the model in simulation

% init comp
[~, ~, n, ~, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

[~, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
labeled_cells = valid;

region_sets = num2cell(valid);
region_area = cell_area(valid);
region_intensity = exp(cell_log_intensity(valid));
num_region = length(region_sets);
region_num_cells = ones(num_region, 1);

[sets_all, log_like_all] = merge_region_fast(num_region, region_area, ...
    region_intensity, region_sets, adj_mat(valid, valid), region_num_cells, n);

% the factor is 4 instead of 6
BIC_all = -2*log_like_all+4*(num_region-1:-1:0)'*log(n);
[~, index_BIC] = min(BIC_all);

selected = sets_all{index_BIC};
n_region = 0;
pred_class_all = zeros(n, 1);
for i = 1:num_region
    if ~isempty(selected{i})
        n_region = n_region+1;
        pred_class_all(selected{i}) = n_region;
    end
end

end
