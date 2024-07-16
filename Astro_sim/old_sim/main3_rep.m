% numerical experiment 3
% several round-shaped extended sources with different intensities

addpath(genpath('/home/minjay/G-SRG'))

rng(1)
T = 500;
n_source = zeros(T, 1);
loc = [0.3 0.3; 0.5 0.7; 0.7 0.3];
radius = [0.15 0.1 0.05];
num_of_photon = num_of_photon_one*ones(1, 3);
metric_all = cell(T, 1);

parfor t = 1:T
    disp(['The ', num2str(t), '-th repetition'])
    % generate simulated data (inhomogeneous Poisson point process)
    X = sim_inhomo_Pois_const([0 1], [0 1], lambda, loc, radius, num_of_photon);

    % init comp
    [cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1]);

    % get seeds
    [seeds, num, invalid] = get_seeds_sim(0.1, 0.9, 0.1, 0.9,...
        0.1, 0.1, 5, cell_log_intensity, cell_area, cx, cy, 2);
    disp(['Number of regions is ', num2str(num)])

    % make a copy of variable seeds
    region_sets = seeds;

    % graph-based SRG
    adj_mat = get_adj_mat( E, n );
    [region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

    [sets_all, log_like_all] = merge_region(num, cell_area, ...
        cell_log_intensity, region_sets, adj_mat, n);

    BIC_all = -2*log_like_all+6*(num-1:-1:0)'*log(n);
    [~, index_BIC] = min(BIC_all);
    n_source(t) = num-index_BIC;
    
    % we assume that the number of sources is correct
    selected = sets_all{num-length(radius)};
    % remove empty cell array elements
    selected = selected(~cellfun(@isempty, selected));
    
    [dr, far, err, res_source_area_sort, source_x, source_y, n_S_correct] =...
        perf_eval(length(radius), loc, radius, cx, cy, selected, cell_log_intensity, cell_area);
    metric_all{t} = [dr far err res_source_area_sort source_x source_y n_S_correct];
end

save(['main3_rep_lambda_', num2str(lambda), '_num_', num2str(num_of_photon_one), '.mat'],...
    'metric_all', 'n_source', 'loc', 'radius', 'num_of_photon')
