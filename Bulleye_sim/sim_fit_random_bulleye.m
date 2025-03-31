
function [labeled_cells, pred_class_all, n_region,rad] = sim_fit_random_bulleye(X, seed, show_plot,gridsize,seedsize,rep_itr)
% fit the model in simulation
rng(seed)
% init comp
rand_num = 3;
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

[region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);

selected_rep = cell(1,rep_itr);
BIC_rep = zeros(1,rep_itr);

%repeat multiple times
parfor ii = 1:rep_itr
    [sets_all, log_like_all] = merge_region_random_fast(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n,rand_num);
    % the factor is 4 instead of 6
    BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
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

%do another beam search merge until two regions left
rad = 0;
selected_nonempty = selected(~cellfun(@isempty,selected));
if n_region==1
    rad = 0;
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
    [seg_contour,~] = get_contour(n, DT,selected_post,adj_mat,invalid,cell_area);
    if length(seg_contour.background)==2
        rad = 0;
    else
        [contourx,contoury] = get_curve(seg_contour.contourV{setdiff(seg_contour.regionId,seg_contour.background)},false);
        [centx,centy] = centroid(polyshape(contourx,contoury));
        [x_sample,y_sample] = sample_curve(contourx,contoury,300,centx, false);
        coef = FD(x_sample,y_sample);
        rad = abs(coef(2,:));
    end
end

if show_plot
    figure
    GRAY = [0.6 0.6 0.6];
    colors = lines(n_region);
    triplot(DT, 'Color', GRAY)
    hold on
    index = 0;
    for i = 1:n_region
        if ~isempty(selected_nonempty{i})
            index = index+1;
            scatter(cx(selected_nonempty{i}), cy(selected_nonempty{i}), 12, colors(index, :), 'filled')
        end
    end
    axis image
    figure
    GRAY = [0.6 0.6 0.6];
    colors = lines(n_region);
    triplot(DT, 'Color', GRAY)
    hold on
    index = 0;
    for i = 1:n_region
        if ~isempty(selected_post{i})
            index = index+1;
            scatter(cx(selected_post{i}), cy(selected_post{i}), 12, colors(index, :), 'filled')
        end
    end
    axis image
end

end
