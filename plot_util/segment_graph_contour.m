function [seg_contour,min_BIC] = segment_graph_contour(X,draw_s,draw_c, Color,N_r,bound,drop,m) 
%
%
%Use SRG-graph to segment the points, merge them based on BIC criterion and
%return the contour of the segmentation. 
%Input: 
% X: all of the points
% draw_s: logic; if true then draw the segmentation
% draw_c: logic; if true then draw the contour
% Color: the color of contour
% N_r: the expected number or bound over number of region  
% bound: N_r is a bound if true, a fixed number if false.
% drop: logic, if true then drop this simulation when the segmentation 
% not equal to true numbe of region
% return segmentation of N_r regions 
% m: panelty term 
%
%Output:
% C: cell array, each element contains x1,x2,y1,y2 representing the 
% contour of of one region 

disp('using the right one')

% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
%adj_mat = get_adj_mat( E, n );

% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[seeds, seeds_rej, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(0.1, 0.9, 0.1, 0.9,...
    0.2, 0.2, 20, cell_log_intensity, cell_area, cx, cy, 2, 20, 5, invalid,adj_mat);
num = num_s+num_s_pt;

% plot the seeds

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
%tic
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

[sets_all, log_like_all] = merge_region(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n);
%toc

BIC_all = -2*log_like_all+m*(num-1:-1:0)'*log(n);
%If there is a restriction, then force to merge to the N_r number of
%regions or within certain bound. 
L = length(BIC_all); 
if isempty(N_r)
    [min_BIC, min_index] = min(BIC_all); 
else
    if ~bound && N_r<=L
        min_index = L+1-N_r;
        min_BIC = BIC_all(min_index);
    elseif bound && N_r<=L
        [min_BIC, min_index] = min(BIC_all(L+1-N_r:end));
        min_index = min_index+L-N_r;
    end
end
disp(min_BIC)
selected = sets_all{min_index};
if draw_s == true
    figure
    %colors = distinguishable_colors(num);
    colors = lines(num);
    plot_segmentation(DT,  min_index, sets_all, cx, cy, colors);
    axis([0 1 0 1])   
    axis square;
    hold on;
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'yticklabel',[])
end

[seg_contour,~] = get_contour(n, DT, selected,adj_mat,invalid,cell_area);


%drop if the segmentation touches the boundary
% if length(unique([x1,x2])) > length(x1) && drop==true
%     x1 = [];
%     x2 = [];
%     y1 = [];
%     y2 = [];
% end
% 
% if N_s ~= N_r && drop == true
%     %disp('the simulation is dropped')
%     x1 = [];
%     x2 = [];
%     y1 = [];
%     y2 = [];
% end




