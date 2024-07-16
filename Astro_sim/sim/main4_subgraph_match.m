% Circular extended source
clear
close all
addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))
addpath(genpath('~/Desktop/gsrg/plot_util'))


GRAY = [0.6 0.6 0.6];
bw = 200;
loc = [0.5 0.5; 0.4 0.4; 0.4 0.6; 0.6 0.4; 0.6 0.6];
radius = [0.25 0.025*ones(1, 4)];

base_num_in_circle = [10 ones(1, 4)];

factor = 30;
sample_factor = 1;
lambda = 1000;
seed = 9;

X = sim_inhomo_Pois_const([0 1], [0 1], lambda, loc, radius, factor * base_num_in_circle, seed);

h = figure;

subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);

%specify the line width of the circles
lw = 2;

subplot(1, 3, 1)

plot_circles(loc, radius,lw)
hold on
scatter(X(:, 1), X(:, 2), 'k.')
axis([0 1 0 1])
axis square
box on
set(gca, 'fontsize', 12)

% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
seedsizes = [1,2,3,4,5];
gridsizes = [0.2,0.15,0.1];
param_num = 0;
for ii = 1:5
    for jj = 1:3
    seedsize = seedsizes(ii);
    gridsize = gridsizes(jj);
    param_num = param_num+1;
    [seeds, seeds_rej, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(0.1, 0.9, 0.1, 0.9,...
        gridsize, gridsize, seedsize, cell_log_intensity, cell_area, cx, cy, 2, 50, seedsize, invalid, adj_mat);
    num = num_s+num_s_pt;
    disp(['Number of regions is ', num2str(num)])
    
    % plot the seeds
    figure
    subplot(1, 3, 2)
    % specify the colormap
    colors = lines(num);
    plot_seeds2(DT, cx, cy, seeds, seeds_pt, seeds_rej, colors, num_s, num_s_pt)
    set(gca, 'fontsize', 12)
    
    seeds_all = [seeds seeds_pt];
    
    % make a copy of variable seeds
    region_sets = seeds_all;
    
    % graph-based SRG
    [region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');
    [region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);

    tic
    [~, ~, sets_all, log_like_all, ~,~, ~] = merge_region_fast_beamsearch_optimization(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n, bw);

    toc
    BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
    [min_BIC, index_BIC] = min(BIC_all);
    subplot(1, 3, 3)
    %specify the color map for background: gray, point: blue, extend: red
    backg_seg = 1;
    extend_seg = 6;
    point_seg = [33,36,37,38];
    colors = lines(num);
    
    plot_segmentation(DT, index_BIC, sets_all, cx, cy, colors)
    hold on
    plot_circles(loc, radius,lw)
    axis image
    set(gca, 'fontsize', 12)
    
    set(h, 'Position', [0, 0, 800, 260]);
    
    selected = sets_all{index_BIC};
    selected_nonempty = selected(~cellfun(@isempty,selected));
    adj_mat_post{param_num} = get_adj_mat_segment(selected_nonempty,n);
    adj_mat_post{param_num} = adj_mat_post{param_num}&adj_mat;
    subplot(1,3,1)
    idx = find(any(adj_mat_post{param_num}));
    scatter(X(idx,1),X(idx,2))
    axis equal
    end
end
uncertainty = zeros(1,n);


for i = 1:param_num
    for j = (i+1):param_num
        for k = 1:n
            
            uncertainty(k) = uncertainty(k) + any(adj_mat_post{i}(k,:)~=adj_mat_post{j}(k,:));
        end
    end
end

figure
pointsize = 20;
scatter(X(:,1), X(:,2), pointsize, uncertainty,'filled');
colormap turbo
hold on
scatter(X(159,1),X(159,2))
find(abs(X(:,1)-0.3956)<0.0001)
adjk = adj_mat_post{i}(159,:);
find(adj_mat(k,:)==1)