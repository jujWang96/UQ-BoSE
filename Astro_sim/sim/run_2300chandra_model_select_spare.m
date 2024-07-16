

close all
addpath(genpath('/Users/joycewang/Desktop/gsrg/NGC2300'))
addpath(genpath('/Users/joycewang/Desktop/gsrg/boot_util'))

load('ngc2300_box_058kev_evt2.mat')
n = length(X);

X = double(unique(X,'rows'));
GRAY = [0.6,0.6,0.6];
P = 14;
threshold = 15;
M = 150;

rand_num = 3;
rep_itr = 5000;
penalty = 6;

seedsize = 5;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[V, R] = voronoiDiagram(DT);

% get seeds
rng(15)%rng(param_idx)
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
%rng(14)
parfor ii = 1:rep_itr

    [sets_all, log_like_all]  = merge_region_random_fast(num_region, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n, rand_num)
    BIC_all = -2*log_like_all+penalty*(num_region-1:-1:0)'*log(n);
    [min_BIC, index_BIC] = min(BIC_all);
    selected = sets_all{index_BIC};

    selected_rep{ii} = selected;
    BIC_rep(ii) = min_BIC;
end
[min_BIC, index_rep] = min(BIC_rep);
selected = selected_rep{index_rep};
selected_nonempty = selected(~cellfun(@isempty,selected));
num_region_nonempty = length(selected_nonempty);
val = get_metric_value_post_seg(n,num_region_nonempty, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
-2*val+penalty*(num_region_nonempty-1)*log(n)
figure
plot_segmentation_wo_voronoi(DT, selected, cx, cy, colors,12,false)

%force merge to two segment
[region_intensity, region_area, region_num_cells, adj_mat_region] =...
        get_region_int_connect(num_region_nonempty, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
[~,~,sets_all_post, log_like_all_post,~,~,~] = merge_region_fast_beamsearch_optimization(num_region_nonempty, region_area, ...
    region_intensity, selected_nonempty, adj_mat_region, region_num_cells, n,30);
selected_post = sets_all_post{num_region_nonempty-1};
figure
plot_segmentation_wo_voronoi(DT, selected_post, cx, cy, colors,12,false)

%iteratively merge the segments and extract the nested boundaries.
selected_nonempty{1} = [selected_nonempty{1},selected_nonempty{7}];
selected_nonempty{7} = [];
figure
plot_segmentation_wo_voronoi(DT, selected_rep{index_rep}, cx, cy, colors,12,false)

structure = [{1},{[1,2]},{[2,3]},{[3,4]},{[4,5]},{5}];
target_x = cell(1,num_region_nonempty-2);
target_y = cell(1,num_region_nonempty-2);

[glb_contour,~] = get_contour(n, DT,selected_nonempty,adj_mat,invalid,cell_area);
contourV = glb_contour.contourV;

[target_x{1},target_y{1}] = get_curve(contourV{5},false);
selected_nonempty{4} = [selected_nonempty{4},selected_nonempty{5}];
selected_nonempty{5} = [];
figure
plot_segmentation_wo_voronoi(DT, selected_nonempty, cx, cy, colors,12,false)

[glb_contour,~] = get_contour(n, DT,selected_nonempty,adj_mat,invalid,cell_area);
contourV = glb_contour.contourV;
[target_x{2},target_y{2}] = get_curve(contourV{4},false);
selected_nonempty{4} = [selected_nonempty{4},selected_nonempty{6}];
selected_nonempty{6} = [];
figure
plot_segmentation_wo_voronoi(DT, selected_nonempty, cx, cy, colors,12,false)

[glb_contour,~] = get_contour(n, DT,selected_nonempty,adj_mat,invalid,cell_area);
contourV = glb_contour.contourV;
[target_x{3},target_y{3}] = get_curve(contourV{4},false);
selected_nonempty{3} = [selected_nonempty{3},selected_nonempty{4}];
selected_nonempty{4} = [];
figure
plot_segmentation_wo_voronoi(DT, selected_nonempty, cx, cy, colors,12,false)

[glb_contour,~] = get_contour(n, DT,selected_nonempty,adj_mat,invalid,cell_area);
contourV = glb_contour.contourV;
[target_x{4},target_y{4}] = get_curve(contourV{3},false);
selected_nonempty{2} = [selected_nonempty{2},selected_nonempty{3}];
selected_nonempty{3} = [];
figure
plot_segmentation_wo_voronoi(DT, selected_nonempty, cx, cy, colors,12,false)

[glb_contour,~] = get_contour(n, DT,selected_nonempty,adj_mat,invalid,cell_area);
contourV = glb_contour.contourV;
[target_x{5},target_y{5}] = get_curve(contourV{2},false);



%select fd pairs with aic/bic criteria
out_bd_x = [0,1,1,0]';
out_bd_y = [0,0,1,1]';
bound = [20,20,20,20,20]; %set upper bound for model selection 


[smooth_polyset,min_pair_num] = model_selection_multiple_fast(X,target_x,target_y,structure,out_bd_x,out_bd_y,["AIC","BIC"],bound);


target_idx = 4;
save('model_select_2300chandra.mat')
