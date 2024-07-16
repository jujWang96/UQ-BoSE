
close all
clear
addpath(genpath('~/src/gsrg/NGC2300XMM/'))
addpath(genpath('~/src/gsrg/NGC2300'))
addpath(genpath('~/src/gsrg/gcr_util'))
addpath(genpath('~/src/gsrg/plot_util'))
addpath(genpath('~/src/gsrg/Astro_sim/util'))
addpath(genpath('~/src/gsrg/G-SRG'))
addpath(genpath('~/src/gsrg/contour_util'))
addpath(genpath('~/src/gsrg/FD_util'))
cd '~/src/gsrg/Astro_sim/real_img/'
load('~/src/gsrg/NGC2300XMM/ngc2300_MOS1_evt_0.5-8.0keV_scaled.mat')

%X = double(unique(X,'rows'));
n = length(X);
% fraction = 1
% id = randsample(n,round(n*fraction));
% X = X(id,:);
X = double(unique(X,'rows'));
GRAY = [0.6,0.6,0.6];
P= 14;
threshold = 20;
M = 150;

rand_num = 5;
rep_itr = 3000;
penalty = 6;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[V, R] = voronoiDiagram(DT);

% get seeds
rng(16) %param_idx
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
delete(gcp);

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
    region_intensity, selected_nonempty, adj_mat_region, region_num_cells, n,20);
selected_post = sets_all_post{num_region_nonempty-1};
figure
plot_segmentation_wo_voronoi(DT, selected_post, cx, cy, colors,12,false)
[seg_contour,~] = get_contour_flexible(n, DT,selected_nonempty,adj_mat,invalid,cell_area);
[target_contour,~] = get_contour_flexible(n, DT,selected_post,adj_mat,invalid,cell_area);
target_contour.contourV
target_id = setdiff(target_contour.regionId,target_contour.background);
if ~length(target_id)==1
    disp("multiple or no target curve identified, plot the segmentation to check")
end
[cxs_t,cys_t] = get_curve_flexible(target_contour.contourV{target_id},false);
target_x  =cxs_t{1}; target_y = cys_t{1};
plot_results = true;
candidates = 2:300  ;
[candidate_aic,candidate_bic,pgon_confine] = model_selection_flexible_fast(X,target_x, target_y, seg_contour,candidates,plot_results);
[centx,centy] = centroid(polyshape(target_x,target_y));
[x_sample,y_sample] = sample_curve(target_x,target_y,300,centx, false);
f_coefs = FD(x_sample,y_sample);

[invx,invy] = iFD(f_coefs,candidate_aic);
figure
plot(pgon_confine)
hold on
plot(invx,invy)
save('flexible_model_select_result_2300XMM.mat','candidate_aic','candidate_bic','pgon_confine')


% pgon_target = polyshape(cxs_t{1},cys_t{1});
% area_target = area(pgon_target);
% x_curves = {};
% y_curves = {};
% pgon_list = {};
% outer_enclose_id = 0;
% outer_match_idx = 0;
% inner_enclose_id = 0;
% inner_match_idx = 0;
% n_sub = 0;
% for objId = seg_contour.regionId
%     [cxs,cys] = get_curve_flexible(seg_contour.contourV{objId},false);
%     pgon_i = polyshape(cxs,cys);
%     max_area_i = 0;
%     pgon_list{end+1} =pgon_i;
%     x_curves{end+1}= cxs;
%     y_curves{end+1} = cys;
%     is_match = false;
%     matched_idx = 0;
%     for i =1:length(cxs)
%         max_area_i = max(max_area_i, polyarea(cxs{i},cys{i}));
%         is_match_i = isequal(sort([cxs{i},cys{i}], 1), sort(pgon_target.Vertices, 1));
%         if is_match_i
%             matched_idx = i;
%         end
%         is_match = is_match || is_match_i;
%         
%         
%     end
%     if is_match && (round(max_area_i,3)> round(area_target,3))
%         outer_enclose_id = objId;
%         outer_match_idx = matched_idx;
%         n_sub = n_sub + sum(isinterior(pgon_i,X(:,1),X(:,2)));
%     end
%     if is_match && (round(max_area_i,3) == round(area_target,3))
%         inner_enclose_id = objId;
%         inner_match_idx = matched_idx;
%         n_sub = n_sub + sum(isinterior(pgon_i,X(:,1),X(:,2)));
%     end
% end
% 
% if (outer_enclose_id==0) || (inner_enclose_id==0)
%     disp("concentric condition is not satisfied")
% end
% 
% %perform model selection 
% [centx,centy] = centroid(polyshape(target_x,target_y));
% [x_sample,y_sample] = sample_curve(target_x,target_y,300,centx, false);
% f_coefs = FD(x_sample,y_sample);
% 
% candidates = 2:30;
% aic_val = NaN(1, length(candidates));
% bic_val = NaN(1, length(candidates));
% 
% for i = 1:length(candidates)
%     candidate = candidates(i);
%     [invx,invy] = iFD(f_coefs,candidate);
%     
%     x_outer_enclose = x_curves{outer_enclose_id};
%     x_outer_enclose{outer_match_idx} = invx;
%     y_outer_enclose = y_curves{outer_enclose_id};
%     y_outer_enclose{outer_match_idx} = invy;
%     
%     pgon_outer = polyshape(x_outer_enclose,y_outer_enclose);
%     
%     x_inner_enclose = x_curves{inner_enclose_id};
%     x_inner_enclose{inner_match_idx} = invx;
%     y_inner_enclose = y_curves{inner_enclose_id};
%     y_inner_enclose{inner_match_idx} = invy;
%     
%     pgon_inner = polyshape(x_inner_enclose,y_inner_enclose);
%     if pgon_outer.NumRegions>1 || pgon_inner.NumRegions>1
%        disp("disgard for intersecting bound") 
%        continue
%     end
%     a_in = area(pgon_inner);
%     n_in = sum(isinterior(pgon_inner,X(:,1),X(:,2)));
%     a_out = area(pgon_outer);
%     n_out = sum(isinterior(pgon_outer,X(:,1),X(:,2)));
%     log_like = n_in*log(n_in/a_in)+n_out*log(n_out/a_out)-sum(log(1:n_sub))-n_sub;
%     aic_val(i) = -2*log_like+(candidate*2+1)*2;
%     bic_val(i) = -2*log_like+(candidate*2+1)*log(n_sub);
% end
% 
% line(candidates, aic_val)
% line(candidates, bic_val)
% 
% for i = 1:length(pgon_list)
%     is_match = false;
%     pgon_i = pgon_list{i};
%     if length(pgon_i.Vertices)==1
%         is_match_j = isequal(sort(pgon_i.Vertices, 1), sort(pgon_target.Vertices, 1));
%     else
%         for j = 1:length(pgon_i.Vertices)
%             is_match_j = isequal(sort(pgon_i.Vertices{j}, 1), sort(pgon_target.Vertices, 1));
%             if is_match_j
%                 plot(pgon_i)
%             end
%         end
%     end
% %     if ~isempty(intersection_shape)
% %         figure
% %         plot(intersection_shape)
% %     end
% end
% 
% pol = polyshape(x_curves,y_curves);
% plot(pol)
% 
% target_x = cxs_t{1};
% target_y = cys_t{1};
%  [centx,centy] = centroid(polyshape(target_x,target_y));
%     [x_sample,y_sample] = sample_curve(target_x,target_y,300,centx, false);
%     f_coefs = FD(x_sample,y_sample);
%     [invx,invy] = iFD(f_coefs,20);
% % %iteratively merge to extract boundaries
% % target_x = {};
% % target_y = {};
% % selected_nonempty_merge = selected_nonempty;
% % seg_countour = get_contour(n, DT,selected_nonempty_merge,adj_mat,invalid,cell_area);
% % 
% % contourV = seg_countour.contourV;
% % [target_x{1},target_y{1}] = get_curve(contourV{1},false);
% % selected_nonempty_merge{1} = [selected_nonempty_merge{1},selected_nonempty_merge{2}];
% % selected_nonempty_merge{2} = [];
% % seg_countour = get_contour(n, DT,selected_nonempty_merge,adj_mat,invalid,cell_area);
% % 
% % contourV = seg_countour.contourV;
% % [target_x{2},target_y{2}] = get_curve(contourV{1},false);
% % selected_nonempty_merge{1} = [selected_nonempty_merge{1},selected_nonempty_merge{3}];
% % 
% % selected_nonempty_merge{3} = [];
% % seg_countour = get_contour(n, DT,selected_nonempty_merge,adj_mat,invalid,cell_area);
% % 
% % contourV = seg_countour.contourV;
% % [target_x{3},target_y{3}] = get_curve(contourV{1},false);
% % 
% % selected_nonempty_merge{1} = [selected_nonempty_merge{1},selected_nonempty_merge{4}];
% % 
% selected_nonempty_merge{4} = [];
% seg_countour = get_contour(n, DT,selected_nonempty_merge,adj_mat,invalid,cell_area);
% 
% contourV = seg_countour.contourV;
% [target_x{4},target_y{4}] = get_curve(contourV{1},false);
% 
% 
% structure = [{1},{[1,2]},{[2,3]},{[3,4]},{4}];
% target_idx = 3;
% out_bd_x = [0,1,1,0]';
% out_bd_y = [0,0,1,1]';
% bd = [20,20,20,20];
% [smooth_polyset,min_pair_num] = model_selection_multiple_fast(X,target_x,target_y,structure,out_bd_x,out_bd_y,["AIC","BIC"],bd);
% 
% [centx,centy] = centroid(polyshape(target_x{target_idx},target_y{target_idx}));
% [x_sample,y_sample] = sample_curve(target_x{target_idx},target_y{target_idx},300,centx, false);
% fourier_coef = FD(x_sample,y_sample)';

