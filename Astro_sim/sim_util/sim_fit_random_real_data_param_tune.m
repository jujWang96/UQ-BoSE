

function [P,threshold] = sim_fit_random_real_data_param_tune(X, seed, show_plot, Ps, thresholds, rand_num, rep_itr, X_obs, min_modelselect,penalty, M)
X = double(unique(X,'rows'));

% fit the model in simulation
rng(seed)
% init comp
BIC_tune = zeros(1,length(Ps)*length(thresholds));
BIC_tune_post = zeros(1,length(Ps)*length(thresholds));

[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

% get seeds
[invalid, ~] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[V, R] = voronoiDiagram(DT);
idx = 0;
for P=Ps
    for threshold = thresholds
        idx = idx+1;
[seeds, num_region] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat,cell_area,P, threshold);
num = num_region;
%disp(['Number of regions is ', num2str(num)])

% make a copy of variable seeds
region_sets = seeds;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

[region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);

if num> M
    [sets_all_greedy, log_like_all_greedy]  = merge_region_fast(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n);
    sets_greedy = sets_all_greedy{num-M+1};
    region_sets = sets_greedy(~cellfun(@isempty,sets_greedy));
    num = length(region_sets); %should be equal to M
    [region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);
end
selected_rep = cell(1,rep_itr);
BIC_rep = zeros(1,rep_itr);

%repeat multiple times
parfor ii = 1:rep_itr
    [sets_all, log_like_all] = merge_region_random_fast(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n,rand_num, ii);
    BIC_all = -2*log_like_all+penalty*(num-1:-1:0)'*log(n);
    [min_BIC, index_BIC] = min(BIC_all);
    selected = sets_all{index_BIC};
    selected_rep{ii} = selected;
    BIC_rep(ii) = min_BIC;
end

[min_BIC, index_rep] = min(BIC_rep);
selected = selected_rep{index_rep};
           

n_region = 0;
pred_class_all = zeros(n, 1);
for i = 1:num
    if ~isempty(selected{i})
        n_region = n_region+1;
        pred_class_all(selected{i}) = n_region;
    end
end

selected_nonempty = selected(~cellfun(@isempty,selected));
num_region_nonempty = length(selected_nonempty);

val = get_metric_value_post_seg(n,num_region_nonempty, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
BIC_tune(idx) = -2*val+penalty*(num_region_nonempty-1)*log(n);
  %do another beam merge until two regions left

if n_region==1
    return 
else
    if n_region==2
        selected_post = selected_nonempty;
    else
    [region_intensity, region_area, region_num_cells, adj_mat_region] =...
        get_region_int_connect(n_region, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
    [~,~,sets_all_post, log_like_all_post,~,~,~] = merge_region_fast_beamsearch_optimization(n_region, region_area, ...
        region_intensity, selected_nonempty, adj_mat_region, region_num_cells, n, 20);
    selected_post = sets_all_post{n_region-1};
    end
    %get the countour and calculate fd
     selected_post_nonempty = selected_post(~cellfun(@isempty,selected_post));
    num_region_post_nonempty = length(selected_post_nonempty);
     val = get_metric_value_post_seg(n,num_region_post_nonempty, cell_area, cell_log_intensity, selected_post_nonempty, adj_mat);
    BIC_tune_post(idx) = -2*val+penalty*(num_region_post_nonempty-1)*log(n);
  
end
    end
end
[~,param_idx] = min(BIC_tune_post);
P = Ps(ceil(param_idx/length(thresholds)));
threshold = thresholds((mod(param_idx-1,length(thresholds)))+1);



end
