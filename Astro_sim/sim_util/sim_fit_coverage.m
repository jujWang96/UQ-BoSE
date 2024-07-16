
function [contourx,contoury,coef] = sim_fit_coverage(X, seed, show_plot,gridsize,seedsize,rep_itr, merge_method)
% fit the model in simulation
if nargin == 6
    merge_method = "random";
end
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
seeds_all = seeds ;

% [seeds, seeds_rej, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(0.1, 0.9, 0.1, 0.9,...
%     0.1, 0.1, 5, cell_log_intensity, cell_area, cx, cy, 2, 50, 5, invalid, adj_mat);
% num = num_s+num_s_pt;
% seeds_all = [seeds seeds_pt];

%disp(['Number of regions is ', num2str(num)])


% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');
colors = lines(length(region_sets));

if show_plot
    fig = figure;
    GRAY = [0.6 0.6 0.6];

    triplot(DT, 'Color', GRAY)
    hold on
    for i = 1:length(region_sets)
        scatter(cx(region_sets{i}), cy(region_sets{i}), 12,  colors(i, :), 'filled')
    end
    
    axis image
end
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

if merge_method == "random"
    selected_rep = cell(1,rep_itr);
    BIC_rep = zeros(1,rep_itr);
    
    %repeat multiple times and select the results with 2 regions 
    parfor ii = 1:rep_itr
        [sets_all, log_like_all] = merge_region_random_fast(num, region_area, ...
        region_intensity, region_sets, adj_mat_region, region_num_cells, n,rand_num, ii );
        % the factor is 4 instead of 6
        BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
        %[min_BIC, index_BIC] = min(BIC_all);
        index_BIC = length(BIC_all)-1; %force merge to 2 regions
        selected = sets_all{index_BIC};
        selected_rep{ii} = selected;
        %BIC_rep(ii) = min_BIC;
        BIC_rep(ii) = BIC_all(index_BIC);
    end
    
    [~, index_rep] = min(BIC_rep);
    disp(BIC_rep(index_rep))
    selected = selected_rep{index_rep};
else
    [sets_all, log_like_all] = merge_region_fast(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n);

    % the factor is 4 instead of 6
    BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
    %[~, index_BIC] = min(BIC_all);
    index_BIC = length(BIC_all)-1;
    disp(BIC_all(index_BIC))

    selected = sets_all{index_BIC};
end
n_region = 0;
pred_class_all = zeros(n, 1);
for i = 1:num
    if ~isempty(selected{i})
        n_region = n_region+1;
        pred_class_all(selected{i}) = n_region;
    end
end

selected_nonempty = selected(~cellfun(@isempty,selected));

selected_post = selected_nonempty;

%     [region_intensity, region_area, region_num_cells, adj_mat_region] =...
%         get_region_int_connect(n_region, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
%     [~,~,sets_all_post, log_like_all_post,~,~,~] = merge_region_fast_beamsearch_optimization(n_region, region_area, ...
%         region_intensity, selected_nonempty, adj_mat_region, region_num_cells, n, 20);
%     selected_post = sets_all_post{n_region-1};

%get the countour and calculate fd
[seg_contour,~] = get_contour(n, DT,selected_post,adj_mat,invalid,cell_area);
if length(seg_contour.background)==2
    coef = zeros(1,300);
    contourx = zeros(1,300);
    contoury = zeros(1,300);
else
    [contourx,contoury] = get_curve(seg_contour.contourV{setdiff(seg_contour.regionId,seg_contour.background)},false);
    [centx,centy] = centroid(polyshape(contourx,contoury));
    [x_sample,y_sample] = sample_curve(contourx,contoury,300,centx, false);
    coef = FD(x_sample,y_sample);

end


if show_plot
    subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);

    figure

    subplot(1,3,1)
    scatter(X(:,1),X(:,2),'k.')
    box on 
    axis equal
    axis([0,1,0,1])
    xlabel(texlabel('X'))
    ylabel(texlabel('Y'))
    set(get(gca,'ylabel'),'rotation',0)
    set(gca,'fontsize',14)
    xticks(0:0.2:1)
    yticks(0:0.2:1)

    subplot(1,3,2)
    %move the background to the beginning
    temp = selected_post(seg_contour.background);
    selected_post(seg_contour.background) = [];
    selected_post = [temp, selected_post];
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
    xlabel(texlabel('X'))
    set(gca,'fontsize',14)
    xticks(0:0.2:1)
    yticks(0:0.2:1)
    axis image
    subplot(1,3,3)
    plot(contourx,contoury,'r')
    box on 
    axis equal
    axis([0,1,0,1])
    xlabel(texlabel('X'))
    set(gca,'fontsize',14)
    xticks(0:0.2:1)
    yticks(0:0.2:1)
end

end
