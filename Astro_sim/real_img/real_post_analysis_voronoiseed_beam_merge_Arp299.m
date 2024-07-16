load('real_full_voronoiseed_random_Arp299.mat')
addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))
addpath(genpath('~/Desktop/gsrg/plot_util'))
addpath(genpath('~/Desktop/gsrg/Astro_sim/sim_util'))
draw_s = false;
draw_c = false;
bound = [0,1;0,1];
rng(1)
n = length(X);
X = double(unique(X,'rows'));

P=25;
bs = 100;
threshold = 10;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[seeds, num] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat,cell_area,P, threshold);
seeds_all = [seeds];

% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');


[log_like_alls, all_seg, sets_all_opt, log_like_all_opt,~,~,~] = merge_region_beam_search(num, cell_area,cell_log_intensity, region_sets, adj_mat, n,100);

BIC_all = -2*log_like_all_opt+6*(num-1:-1:0)'*log(n);
[min_BIC_bs, index_BIC] = min(BIC_all);
%
selected = sets_all_opt{index_BIC};

% plot the data
fig = figure;
scatter(X(:, 1), X(:, 2), 1,'k.')
axis image
box on
set(gca, 'fontsize', 18)
min_white_margin(gca);
set(fig, 'Position', [0, 0, 490, 380]);

% get seeds
[seeds, num] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat,cell_area,P, threshold);
fig = figure;
plot_seeds(DT, cx, cy, seeds, [], [], lines(num), num, 0)


fig = figure;
triplot(DT, 'Color', GRAY)
hold on
% the final resul

selected_nonempty = {};
index = 0;
for i = 1:length(selected)
    if ~isempty(selected{i})
        index = index + 1;
        selected_nonempty{index} = selected{i};
    end
end

for i = 1:length(selected_nonempty)-1
    for j = i+1:length(selected_nonempty)
        % sorted by area from the largest to the smallest
        if sum(cell_area(selected_nonempty{i})) < sum(cell_area(selected_nonempty{j}))
            tmp = selected_nonempty{i};
            selected_nonempty{i} = selected_nonempty{j};
            selected_nonempty{j} = tmp;
        end
    end
end

for i = 1:length(selected_nonempty)
    log_int = log(sum(exp(cell_log_intensity(selected_nonempty{i})).*cell_area(selected_nonempty{i}))/sum(cell_area(selected_nonempty{i})));
    scatter(cx(selected_nonempty{i}), cy(selected_nonempty{i}), 12, log_int*ones(length(selected_nonempty{i}), 1), 'filled')
end
colorbar('EastOutside')
colormap(hsv)
caxis([10 18])
axis image
set(gca, 'fontsize', 12)
min_white_margin(gca);
%saveas(fig, 'segment_result_points', 'png')
% point sources (to be compared with the result from the wavdetect
% algorithm)
% not to be fancy, choose a cutoff point
cutoff_point_source = 0.0003;
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

fig = figure;
triplot(DT, 'Color', GRAY)
hold on
for i = 1:num_nonempty
    % if area is less than the threshold
    if area_all(i) < cutoff_point_source
        scatter(cx(selected_nonempty{i}), cy(selected_nonempty{i}), 12, log_int_all(i)*ones(length(selected_nonempty{i}), 1), 'filled')
    end
end
colorbar('EastOutside')
colormap(hsv)
axis image
set(gca, 'fontsize', 14)
min_white_margin(gca);


%[glb_contour,glb_region_intensity,glb_selected,glb_min_BIC,seg_contour,seg_region_intensity,seg_selected,seg_min_BIC] = gSRG_random3_voronoiseeds(X,true,true, 'k',[],true,true,2,rand_step_num,seg_rep,full_rep,10,P,threshold);
save('real_full_voronoiseed_beam_Arp299_penalty6.mat','X','selected')