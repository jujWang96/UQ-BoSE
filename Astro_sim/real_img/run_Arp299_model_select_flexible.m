
close all
clear
addpath(genpath('/Users/joycewang/src/gsrg/Arp299XMM'))
addpath(genpath('/Users/joycewang/src/gsrg/contour_util'))
addpath(genpath('/Users/joycewang/src/gsrg/boot_util'))
addpath(genpath('/Users/joycewang/src/gsrg/FD_util'))
addpath '/Users/joycewang/src/gsrg/Astro_sim/real_img'
cd('~/src/gsrg/Arp299XMM/')
load('Arp299_MOS1_evt_0.5-8.0keV_scaled.mat')
load('tune_param_Arp299.mat')
%X = double(unique(X,'rows'));
n = length(X);
% fraction = 1
% id = randsample(n,round(n*fraction));
% X = X(id,:);
X = double(unique(X,'rows'));
GRAY = [0.6,0.6,0.6];
%P= 17;
%threshold = 5;
%M = 150;

rand_num = 3;
rep_itr = 10000;
penalty = 6;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[V, R] = voronoiDiagram(DT);

% get seeds
rng(param_idx) %param_idx
[seeds, num_region] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat,cell_area,P, threshold);

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
sets_greedy = region_sets;

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
selected_nonempty = selected(~cellfun(@isempty,selected));
num_region_nonempty = length(selected_nonempty);
val = get_metric_value_post_seg(n,num_region_nonempty, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
-2*val+penalty*(num_region_nonempty-1)*log(n)
figure
plot_segmentation_wo_voronoi(DT, selected_rep{index_rep}, cx, cy, colors,12,false)

%force merge to two segment
[region_intensity, region_area, region_num_cells, adj_mat_region] =...
        get_region_int_connect(num_region_nonempty, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
[~,~,sets_all_post, log_like_all_post,~,~,~] = merge_region_fast_beamsearch_optimization(num_region_nonempty, region_area, ...
    region_intensity, selected_nonempty, adj_mat_region, region_num_cells, n,100);
selected_post = sets_all_post{num_region_nonempty-1};
figure
plot_segmentation_wo_voronoi(DT, selected_post, cx, cy, colors,12,false)
[seg_contour,~] = get_contour_flexible(n, DT,selected_nonempty,adj_mat,invalid,cell_area);
[target_contour,~] = get_contour_flexible(n, DT,selected_post,adj_mat,invalid,cell_area);
target_contour.contourV
target_id = setdiff(target_contour.regionId,target_contour.background);
for i = seg_contour.regionId
    disp(seg_contour.segArea(i)*seg_contour.segIntensity(i))
end
if ~length(target_id)==1
    disp("multiple or no target curve identified, plot the segmentation to check")
end
[cxs_t,cys_t] = get_curve_flexible(target_contour.contourV{target_id},false);
target_x  =cxs_t{1}; target_y = cys_t{1};
plot_results = true;
candidates = 2:300 ;
[candidate_aic,candidate_bic,pgon_confine] = model_selection_flexible_fast(X,target_x, target_y, seg_contour,candidates,plot_results);
[candidate_aic,candidate_bic,pgon_confine] = model_selection_flexible_robust(X,target_x, target_y, seg_contour,candidates,plot_results);

[centx,centy] = centroid(polyshape(target_x,target_y));
[x_sample,y_sample] = sample_curve(target_x,target_y,300,centx, false);
f_coefs = FD(x_sample,y_sample);

[invx,invy] = iFD(f_coefs,candidate_aic);
figure
plot(pgon_confine)
hold on
plot(invx,invy)

%save('flexible_model_select_result_299XMM.mat','candidate_aic','candidate_bic','pgon_confine')
