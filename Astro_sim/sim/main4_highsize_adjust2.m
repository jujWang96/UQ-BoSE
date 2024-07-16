% Circular extended source
%adjust for the area also but not just the intensity 
clear
close all
addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))
addpath(genpath('~/Desktop/gsrg/plot_util'))
filename = 'main4_lambda2000_';
penalty = 3;
filename = strcat(filename, 'penalty',num2str(penalty),'_');
seed = 1;
filename = strcat(filename, 'seed',num2str(seed));

rep_itr = 1000;
rand_num = 5;

GRAY = [0.6 0.6 0.6];

loc = [0.5 0.5; 0.4 0.4; 0.4 0.6; 0.6 0.4; 0.6 0.6];
radius = [0.25 0.025*ones(1, 4)];

base_num_in_circle = [10 ones(1, 4)];

factor = 30;
sample_factor = 1;
lambda = 2000;
instrument = [0.1,0;1,0.9;0.9,1;0,0.1];
% figure
% plot(polyshape(instrument))
% axis equal
instrument_factor = 1/2;
X = sim_inhomo_Pois_const_adjust([0 1], [0 1], lambda, loc, radius, factor * base_num_in_circle, instrument_factor,instrument,seed,false);


subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);

%specify the line width of the circles
lw = 2;
h = figure;
subplot(1, 3, 1)

% plot_circles(loc, radius,lw)
% hold on
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
[seeds, seeds_rej, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(0.1, 0.9, 0.1, 0.9,...
    0.2, 0.2, 5, cell_log_intensity, cell_area, cx, cy, 2, 50, 5, invalid, adj_mat);
num = num_s+num_s_pt;
disp(['Number of regions is ', num2str(num)])

% plot the seeds
%subplot(1, 3, 2)
% specify the colormap
colors = lines(num);
% plot_seeds2(DT, cx, cy, seeds, seeds_pt, seeds_rej, colors, num_s, num_s_pt)
% set(gca, 'fontsize', 12)

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');
%plot_segmentation_wo_voronoi(DT, region_sets, cx, cy, colors,10,true)

[region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);


rng(seed)
selected_rep = cell(1,rep_itr);
    BIC_rep = zeros(1,rand_num);
    tic;
    for ii = 1:rep_itr
    
        [sets_all, log_like_all]  = merge_region_random_fast(num, region_area, ...
        region_intensity, region_sets, adj_mat_region, region_num_cells, n, rand_num);
        BIC_all = -2*log_like_all+penalty*(num-1:-1:0)'*log(n);
        [min_BIC, index_BIC] = min(BIC_all);
        %index_BIC = length(BIC_all)-1;
        %min_BIC = BIC_all(index_BIC)
        selected = sets_all{index_BIC};
        
        selected_rep{ii} = selected;
        BIC_rep(ii) = min_BIC;
       
    end
    toc;
    [min_BIC_rep, index_rep] = min(BIC_rep);
    selected = selected_rep{index_rep};

subplot(1, 3, 2)
%specify the color map for background: gray, point: blue, extend: red
backg_seg = 1;
extend_seg = 6;
point_seg = [33,36,37,38];
colors = lines(num);%[GRAY;[0.8500,0.3250,0.0980];repmat([0,0.4470,0.7410],10,1)];

plot_segmentation_wo_voronoi(DT, selected, cx, cy, colors, 10,true)
hold on
plot_circles(loc, radius,lw)



%% use adjusted intensity 
% make a copy of variable seeds

adjust_indc = inpolygon(X(:,1),X(:,2),instrument(:,1),instrument(:,2));
%adjust for the number of neighbours that is inside the region 
nbr_in = adj_mat*adjust_indc;
nbr_total = sum(adj_mat)';
nbr_prop = nbr_in./nbr_total;
nbr_prop = ones(size(X,1),1);
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp_adjust(X, [0 1], [0 1],nbr_prop*instrument_factor, adjust_indc,ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[seeds, seeds_rej, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(0.1, 0.9, 0.1, 0.9,...
    0.2, 0.2, 5, cell_log_intensity, cell_area, cx, cy, 2, 50, 5, invalid, adj_mat);
num = num_s+num_s_pt;
disp(['Number of regions is ', num2str(num)])

% plot the seeds
% specify the colormap
colors = lines(num);
% figure
% plot_seeds2(DT, cx, cy, seeds, seeds_pt, seeds_rej, colors, num_s, num_s_pt)
% set(gca, 'fontsize', 12)

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

%cell_log_intensity = cell_log_intensity-log(instrument_factor)*adjust_indc.*nbr_prop;

% graph-based SRG

[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');
% figure
% plot_segmentation_wo_voronoi(DT, region_sets, cx, cy, colors,10,true)
[region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);
%plot intensity 
% figure
% triplot(DT, 'Color', GRAY)
% hold on
% for i = 1:length(region_sets)
%     log_int = log(region_intensity(i));
%     scatter(cx(region_sets{i}), cy(region_sets{i}), 12, log_int*ones(length(region_sets{i}), 1), 'filled')
% end
% axis equal
% axis([0 1 0 1])
% tic
%  [sets_all, log_like_all]  = merge_region_fast_with_step(cx,cy,num, region_area, ...
%     region_intensity, region_sets, adj_mat_region, region_num_cells, n);
% toc
% BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
% [min_BIC, index_BIC] = min(BIC_all);
%  figure
%    plot_segmentation_wo_voronoi(DT, sets_all{index_BIC}, cx, cy, colors,10,true)
rng(seed)

selected_rep = cell(1,rep_itr);
    BIC_rep = zeros(1,rand_num);
    tic;
    for ii = 1:rep_itr
    
        [sets_all, log_like_all]  = merge_region_random_fast(num, region_area, ...
        region_intensity, region_sets, adj_mat_region, region_num_cells, n, rand_num);
        BIC_all = -2*log_like_all+penalty*(num-1:-1:0)'*log(n);
        [min_BIC, index_BIC] = min(BIC_all);
        %index_BIC = length(BIC_all)-1;
        %min_BIC = BIC_all(index_BIC)
        selected = sets_all{index_BIC};
%         figure
%          plot_segmentation_wo_voronoi(DT, selected, cx, cy, colors,10,true)
%         text(0.5,0.5,num2str(min_BIC))
%         text(0.4,0.5,num2str(ii))
        selected_rep{ii} = selected;
        BIC_rep(ii) = min_BIC;
       
    end
    toc;
    [min_BIC_rep, index_rep] = min(BIC_rep);
    selected = selected_rep{index_rep};
    %%check the calculation of BIC is correct 
selected_nonempty = selected(~cellfun(@isempty,selected));
    num_region_nonempty = length(selected_nonempty);
    val = get_metric_value_post_seg(n,num_region_nonempty, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
    BIC_check  =  -2*val+penalty*(num_region_nonempty-1)*log(n);



% figure
% %specify the color map for background: gray, point: blue, extend: red
% backg_seg = 1;
% extend_seg = 6;
% point_seg = [33,36,37,38];
% colors = lines(num);%[GRAY;[0.8500,0.3250,0.0980];repmat([0,0.4470,0.7410],10,1)];
% 
subplot(1,3,3)
 plot_segmentation_wo_voronoi(DT, selected, cx, cy, colors,10,true)
%   text(0.5,0.5,num2str(min_BIC_rep))
%    text(0.4,0.5,"best")
hold on
plot_circles(loc, radius,lw)
axis image
set(gca, 'fontsize', 12)

set(h, 'Position', [0, 0, 800, 260]);
saveas(h,filename,'epsc')