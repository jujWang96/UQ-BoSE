% Circular extended source
clear
close all
addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))
addpath(genpath('~/Desktop/gsrg/plot_util'))
addpath(genpath('~/Desktop/gsrg/contour_util'))

fileID = fopen('arc_10sim_center.txt','w');
formatSpec = '(%4.5f,%4.5f) ';
GRAY = [0.6 0.6 0.6];
colr = [0.5, 0.5, 0];	
pt_sz = 8;
lw_circle = 2;
lw_arc = 2;
loc_ring = [0.3 0.5];
radius_in = 0.2;
radius_out = 0.4;
loc = [0.5 0.7; 0.6 0.5; 0.5 0.3];
radius = 0.025*ones(1, 3);
base_num_in_ring = 10;
base_num_in_circle = ones(1, 3);

factor = 30;
sample_factor = 1;
lambda = 1000;
seed = 0;
rep = 10; %repeat 9 times for uncertainty 
X = sim_inhomo_Pois_const_ring(loc_ring, radius_out, radius_in, factor * base_num_in_ring, seed);
X = [X; sim_inhomo_Pois_const([0 1], [0 1], lambda, loc, radius, factor * base_num_in_circle)];

h = figure;

subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);

%specify the line width of the circles
lw = 2;
lw_seg = 2;

subplot(1, 3, 1)
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
subplot(1, 3, 2)
% specify the colormap
colors = lines(num);
plot_seeds2(DT, cx, cy, seeds, seeds_pt, seeds_rej, colors, num_s, num_s_pt)
hold on
set(gca, 'fontsize', 12)
plot_circles(loc, radius,lw_circle,colr)
plot_arc(loc_ring, radius_in, radius_out,lw_arc,colr)

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);

BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);

subplot(1, 3, 3)
%specify the color map for background: gray, point: blue, extend: red
backg_seg = 1;
extend_seg = 6;
point_seg = [33,36,37,38];
colors = [[0,0,0];[0.8500,0.3250,0.0980];repmat([0,0.4470,0.7410],4,1)];

%plot_segmentation(DT, index_BIC, sets_all, cx, cy, colors)
selected = sets_all{index_BIC};

[area_seg,flux_seg,ratio_seg,x_seq,y_seq] = get_area_flux(X,selected,cell_area,true); 
for i = 1:length(x_seq)
    fprintf(fileID,formatSpec,x_seq(i),y_seq(i));
end
fprintf(fileID,'\n');
% writematrix(ratio_seg,'arc_10sim.txt','WriteMode','append')
% writematrix(area_seg,'arc_10sim.txt','WriteMode','append')
% writematrix(flux_seg,'arc_10sim.txt','WriteMode','append')

DT_pt = DT;
selected_pt = selected;
cx_pt = cx;
cy_pt = cy;

% selected = sets_all{index_BIC};
%     
%     % drop empty entries in the cell array
%     num = length(selected);
%     
%     area_all = [];
%     log_int_all = [];
%     selected_nonempty = {};
%     index = 0;
%     for i = 1:num
%         if ~isempty(selected{i})
%             index = index+1;
%             selected_nonempty{index} = selected{i};
%             area = sum(cell_area(selected{i}));
%             area_all = [area_all, area];
%             log_int = log(sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i}))/area);
%             log_int_all = [log_int_all, log_int];
%         end
%     end
%     
%     num_nonempty = length(selected_nonempty);
%     
%     % get source boundaries formed by Voronoi cell edges
%     [V, R] = voronoiDiagram(DT);
%     vx_edges_all = {};
%     vy_edges_all = {};
%     for j = 1:num_nonempty
%         cells_in_region = R(selected_nonempty{j});
%         edges = [];
%         for i = 1:length(cells_in_region)
%             % each row records two vertex indices of an edge
%             new_edges = [cells_in_region{i}' [cells_in_region{i}(2:end)'; cells_in_region{i}(1)]];
%             % sort vertex indices to avoid ambiguity
%             new_edges = sort(new_edges, 2);
%             edges = [edges; new_edges];
%         end
%         [unique_edges, ~ , ind] = unique(edges, 'rows');
%         % get the count of each unique edge
%         counts = histc(ind, unique(ind));
%         % remove edges that appear at least twice to get the boundary edges
%         edges = unique_edges(counts==1, :);
%         % x-axis
%         vx_edges = [V(edges(:, 1), 1)'; V(edges(:, 2), 1)'];
%         % y-axis
%         vy_edges = [V(edges(:, 1), 2)'; V(edges(:, 2), 2)'];
%         nume = size(vx_edges, 2);
%         vx_edges = [vx_edges; NaN(1, nume)];
%         vx_edges = vx_edges(:);
%         vx_edges_all{j} = vx_edges;
%         vy_edges = [vy_edges; NaN(1, nume)];
%         vy_edges = vy_edges(:);
%         vy_edges_all{j} = vy_edges;
%     end
% 
% vx_edges_all_pt = vx_edges_all;
% 
% vy_edges_all_pt = vy_edges_all;
% 
% num_nonempty_pt = num_nonempty;

%add uncertainty curve 
for ii = 1:rep
    
    seed = ii;
    X = sim_inhomo_Pois_const_ring(loc_ring, radius_out, radius_in, factor * base_num_in_ring, seed);
    X = [X; sim_inhomo_Pois_const([0 1], [0 1], lambda, loc, radius, factor * base_num_in_circle)];

    % init comp
    [cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
    adj_mat = get_adj_mat( E, n );
    
    % get seeds
    [invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
    [seeds, seeds_rej, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(0.1, 0.9, 0.1, 0.9,...
        0.2, 0.2, 5, cell_log_intensity, cell_area, cx, cy, 2, 50, 5, invalid, adj_mat);
    num = num_s+num_s_pt;
    disp(['Number of regions is ', num2str(num)])
    
    seeds_all = [seeds seeds_pt];
    
    % make a copy of variable seeds
    region_sets = seeds_all;
    
    % graph-based SRG
    [region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');
    
    [sets_all, log_like_all] = merge_region(num, cell_area, ...
        cell_log_intensity, region_sets, adj_mat, n);
    
    BIC_all = -2*log_like_all+4*(num-1:-1:0)'*log(n);
    [min_BIC, index_BIC] = min(BIC_all);
    
    
    selected = sets_all{index_BIC};
    [area_seg,flux_seg,ratio_seg,x_seq,y_seq] = get_area_flux(X,selected,cell_area,true); 
    for i = 1:length(x_seq)
        fprintf(fileID,formatSpec,x_seq(i),y_seq(i));
    end
    fprintf(fileID,'\n');%     writematrix(ratio_seg,'arc_10sim.txt','WriteMode','append')
%     writematrix(area_seg,'arc_10sim.txt','WriteMode','append')
%     writematrix(flux_seg,'arc_10sim.txt','WriteMode','append')
% 
%     
    % drop empty entries in the cell array
    num = length(selected);
    
    area_all = [];
    log_int_all = [];
    selected_nonempty = {};
    index = 0;
    for i = 1:num
        if ~isempty(selected{i})
            index = index+1;
            selected_nonempty{index} = selected{i};
            area = sum(cell_area(selected{i}));
            area_all = [area_all, area];
            log_int = log(sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i}))/area);
            log_int_all = [log_int_all, log_int];
        end
    end
    
    num_nonempty = length(selected_nonempty);
    
    % get source boundaries formed by Voronoi cell edges
    [V, R] = voronoiDiagram(DT);
    vx_edges_all = {};
    vy_edges_all = {};
    area_seg = [];
    flux_seg = [];
    for j = 1:num_nonempty
        cells_in_region = R(selected_nonempty{j});
        edges = [];
        for i = 1:length(cells_in_region)
            % each row records two vertex indices of an edge
            new_edges = [cells_in_region{i}' [cells_in_region{i}(2:end)'; cells_in_region{i}(1)]];
            % sort vertex indices to avoid ambiguity
            new_edges = sort(new_edges, 2);
            edges = [edges; new_edges];
        end
        [unique_edges, ~ , ind] = unique(edges, 'rows');
        % get the count of each unique edge
        counts = histc(ind, unique(ind));
        % remove edges that appear at least twice to get the boundary edges
        edges = unique_edges(counts==1, :);
        % x-axis
        vx_edges = [V(edges(:, 1), 1)'; V(edges(:, 2), 1)'];
        % y-axis
        vy_edges = [V(edges(:, 1), 2)'; V(edges(:, 2), 2)'];
        nume = size(vx_edges, 2);
        vx_edges = [vx_edges; NaN(1, nume)];
        vx_edges = vx_edges(:);
        vx_edges_all{j} = vx_edges;
        vy_edges = [vy_edges; NaN(1, nume)];
        vy_edges = vy_edges(:);
        vy_edges_all{j} = vy_edges;
    end
    
    hold on
    for j = 2:num_nonempty
        % original red is brighter
        if ii==8
            continue
        end
        line(vx_edges_all{j}, vy_edges_all{j}, 'Color', GRAY,'LineWidth',lw_seg)
    end
end

plot_segmentation_wo_voronoi(DT_pt, selected_pt, cx_pt, cy_pt, colors,pt_sz,false) %default pt_size = 12
box on

% for j = 2
%     % original red is brighter
%     line(vx_edges_all_pt{j}, vy_edges_all_pt{j}, 'Color', [0.8500,0.3250,0.0980],'LineWidth',lw_seg)
% end
% for j = 3:num_nonempty_pt
%     % original red is brighter
%     line(vx_edges_all_pt{j}, vy_edges_all_pt{j}, 'Color', [0,0.4470,0.7410],'LineWidth',lw_seg)
% end
%plot_circles(loc, radius,lw)
axis image
set(gca, 'fontsize', 12)

set(h, 'Position', [0, 0, 800, 260]);
%ct = get_contour(n, DT,sets_all{index_BIC},adj_mat,invalid,cell_area);

%line(ct.contourV{1},ct.contourV{1}(:,3))

%get_curve(ct.contourV{6},true)
