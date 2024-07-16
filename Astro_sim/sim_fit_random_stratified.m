function [contourx,contoury, fourier_coef, raw_contour,selected_nonempty,min_BIC,min_BIC_post] = sim_fit_random_stratified(X, seed, show_plot, P, threshold, rand_num, rep_itr,penalty, M)
% use a stratified sampling seed intialization
% use a randomized merger to merge the oversegmentation with 
% M,rand_num, rep_itr specify the number of random merge step, the
% candidate merger size and the repetition number
%X = unique(X,'rows');

K = 300;
fourier_coef = zeros(1,K);
contourx = zeros(1,K);
contoury = zeros(1,K);

% fit the model in simulation
rng(seed)
% init comp

[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

% get seeds
[invalid, ~] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[seeds, num_region] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat,cell_area,P, threshold);
num = num_region;
%disp(['Number of regions is ', num2str(num)])

% make a copy of variable seeds
region_sets = seeds;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

[region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);

% perform greedy merge first 
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
delete(gcp)

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


%do another beam merge until two regions left
selected_nonempty = selected(~cellfun(@isempty,selected));
[raw_contour,~] = get_contour_flexible(n, DT,selected_nonempty,adj_mat,invalid,cell_area);

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
        BIC_all_post = -2*log_like_all_post+penalty*(n_region-1:-1:0)'*log(n);
        min_BIC_post = BIC_all_post(n_region-1);
        selected_post_nonempty = selected_post(~cellfun(@isempty,selected_post));
        num_region_post_nonempty = length(selected_post_nonempty);
        val = get_metric_value_post_seg(n,num_region_post_nonempty, cell_area, cell_log_intensity, selected_post_nonempty, adj_mat);
        BIC_val = -2*val+penalty*(num_region_post_nonempty-1)*log(n);
        disp(strcat("BIC_post_min:",num2str(min_BIC_post),";BIC_val:",num2str(BIC_val)))
     
    end
    %get the countour and calculate fd
    [seg_contour,~] = get_contour_flexible(n, DT,selected_post,adj_mat,invalid,cell_area);
    
    % could not extract boundary because the central segment touches the
    % boundary and is not complete 
    if length(seg_contour.background)==2

        return
    end
    target_id = setdiff(seg_contour.regionId,seg_contour.background);
    [contourx,contoury] = get_curve(seg_contour.contourV{target_id},false);
    [centx,centy] = centroid(polyshape(contourx,contoury));
    [x_sample,y_sample] = sample_curve(contourx,contoury,300,centx, false);
    fourier_coef = FD(x_sample,y_sample);

end

if show_plot
         GRAY = [0.6 0.6 0.6];

%      figure
%      colors = lines(n_region);
%     triplot(DT, 'Color', GRAY)
%     hold on
%     index = 0;
%     for i = 1:n_region
%         if ~isempty(selected_nonempty{i})
%             index = index+1;
%             scatter(cx(selected_nonempty{i}), cy(selected_nonempty{i}), 12, colors(index, :), 'filled')
%         end
%     end
%     axis image
%     figure
%     GRAY = [0.6 0.6 0.6];
%     colors = lines(n_region);
%     triplot(DT, 'Color', GRAY)
%     hold on
%     index = 0;
%     for i = 1:n_region
%         if ~isempty(selected_post{i})
%             index = index+1;
%             scatter(cx(selected_post{i}), cy(selected_post{i}), 12, colors(index, :), 'filled')
%         end
%     end
%     axis image
    figure 
    hold on
    scatter(X(:,1),X(:,2),2,GRAY,'filled')
    axis equal
    axis([0,1,0,1])
    plot(x_sample,y_sample,'Color','r','LineWidth',2)

end

end
