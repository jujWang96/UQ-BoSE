
close all
addpath(genpath('/Users/joycewang/Desktop/gsrg/NGC2300XMM'))
load('ngc2300_MOS1_evt_0.5-8.0keV_scaled.mat')

%X = double(unique(X,'rows'));
n = length(X);
% fraction = 1
% id = randsample(n,round(n*fraction));
% X = X(id,:);
X = double(unique(X,'rows'));
GRAY = [0.6,0.6,0.6];
Ps=[5,8,11,14,17,20]
thresholds = [5,10,15,20]
factor = 2;
M = 150;
BIC_tune = zeros(1,length(Ps)*length(thresholds));
BIC_tune_group = zeros(1,length(Ps)*length(thresholds));

rand_num = 5;
rep_itr = 3000;
penalty = 12;
seedsize = 5;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[V, R] = voronoiDiagram(DT);
idx = 0;
for P = Ps
    for threshold = thresholds
    % get seeds
    idx = idx+1;
    rng(idx)
    [invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
    [seeds, num_region] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat,cell_area,P, threshold);
    %[seeds, ~,  num_s] = get_seeds_sim(0.1, 0.9, 0.1, 0.9,...
     %      gridsize,gridsize, seedsize, cell_log_intensity, cell_area, cx, cy, 2, invalid,adj_mat);
     %num_region = num_s;
    %disp(['Number of regions is ', num2str(num)])
    
    seeds_all = seeds ;
    
    % make a copy of variable seeds
    
    disp(strcat('The number of seeds is',num2str(num_region)))
    
    
    colors = lines(num_region);
    figure
    triplot(DT, 'Color', GRAY)
    hold
    for i = 1:num_region
        scatter(cx(seeds{i}), cy(seeds{i}), 12, colors(i, :), 's', 'filled')
    end
    axis equal
    
    region_sets = seeds;
    
    % graph-based SRG
    [region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');
    
    [region_intensity, region_area, region_num_cells, adj_mat_region] =...
        get_region_int_connect(num_region, cell_area, cell_log_intensity, region_sets, adj_mat);
    
    %greedy merge until M region left then do randomized search
    if num_region> M
        [sets_all_greedy, log_like_all_greedy]  = merge_region_fast(num_region, region_area, ...
        region_intensity, region_sets, adj_mat_region, region_num_cells, n);
        sets_greedy = sets_all_greedy{num_region-M+1};
        region_sets = sets_greedy(~cellfun(@isempty,sets_greedy));
        num_region = length(region_sets); %should be equal to M
        [region_intensity, region_area, region_num_cells, adj_mat_region] =...
        get_region_int_connect(num_region, cell_area, cell_log_intensity, region_sets, adj_mat);
    end
    selected_rep = cell(1,rep_itr);
    BIC_rep = zeros(1,rand_num);
    tic;
    parfor ii = 1:rep_itr
    
        [sets_all, log_like_all]  = merge_region_random_fast(num_region, region_area, ...
        region_intensity, region_sets, adj_mat_region, region_num_cells, n, rand_num)
        BIC_all = -2*log_like_all+penalty*(num_region-1:-1:0)'*log(n);
        [min_BIC, index_BIC] = min(BIC_all);
        selected = sets_all{index_BIC};
    
        selected_rep{ii} = selected;
        BIC_rep(ii) = min_BIC;
    end
    toc;
    [min_BIC, index_rep] = min(BIC_rep);
    selected = selected_rep{index_rep};
    selected_nonempty = selected(find(~cellfun(@isempty,selected)));
    num_region_nonempty = length(selected_nonempty);
    val = get_metric_value_post_seg(n,num_region_nonempty, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
    BIC_tune(idx) = -2*val+penalty*(num_region_nonempty-1)*log(n);
    
    
    figure
    plot_segmentation(DT, index_rep, selected_rep, cx, cy, colors)
    
    [seeds_group, seeds_rej, num_region] = get_seeds_sim_voronoi_group(cx, cy,invalid, adj_mat,cell_area,P, threshold,seedsize,factor);
    colors = lines(num_region);
    
    figure
    triplot(DT, 'Color', GRAY)
    hold
    for i = 1:num_region
        scatter(cx(seeds_group{i}), cy(seeds_group{i}), 12, colors(i, :), 's', 'filled')
    end
    axis equal
    
    region_sets = seeds_group;
    
    % graph-based SRG
    [region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');
    %figure
    %plot_segmentation(DT, 1, {region_sets}, cx, cy, colors)
    
    [region_intensity, region_area, region_num_cells, adj_mat_region] =...
        get_region_int_connect(num_region, cell_area, cell_log_intensity, region_sets, adj_mat);
    
     %greedy merge until M region left then do randomized search
    if num_region> M
        [sets_all_greedy, log_like_all_greedy]  = merge_region_fast(num_region, region_area, ...
        region_intensity, region_sets, adj_mat_region, region_num_cells, n);
        sets_greedy = sets_all_greedy{num_region-M+1};
        region_sets = sets_greedy(~cellfun(@isempty,sets_greedy));
        num_region = length(region_sets); %should be equal to M
        [region_intensity, region_area, region_num_cells, adj_mat_region] =...
        get_region_int_connect(num_region, cell_area, cell_log_intensity, region_sets, adj_mat);
    end
    
    selected_rep = cell(1,rep_itr);
    BIC_rep = zeros(1,rand_num);
    tic;
    parfor ii = 1:rep_itr
    
        [sets_all, log_like_all]  = merge_region_random_fast(num_region, region_area, ...
        region_intensity, region_sets, adj_mat_region, region_num_cells, n, rand_num)
        BIC_all = -2*log_like_all+penalty*(num_region-1:-1:0)'*log(n);
        [min_BIC, index_BIC] = min(BIC_all);
        selected = sets_all{index_BIC};
    
        selected_rep{ii} = selected;
        BIC_rep(ii) = min_BIC;
    end
    toc;
    [min_BIC, index_rep] = min(BIC_rep);
    selected = selected_rep{index_rep};
    figure
    plot_segmentation(DT, index_rep, selected_rep, cx, cy, colors)
    
    selected_nonempty = selected(~cellfun(@isempty,selected));
    num_region_nonempty = length(selected_nonempty);
    val = get_metric_value_post_seg(n,num_region_nonempty, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
    BIC_tune_group(idx) = -2*val+penalty*(num_region_nonempty-1)*log(n);
    end
end
%the best parameter with single seed picked is P=14,S = 10 (-8.3976)
%the best parameter with group seed picked is P=20,S = 10 (-8.4456)
