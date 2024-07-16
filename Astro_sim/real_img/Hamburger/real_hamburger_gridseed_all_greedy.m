close all
clear
cd('~/Desktop/gsrg/Astro_sim/real_img/Hamburger')
addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))
fid=fopen('Hamburger_11007/hamburger_0123.txt');
C=textscan(fid,'%f%f%f%f%f%f%f%f%f','HeaderLines',11,'delimiter','\r\n');
fclose(fid);

X = [C{1},C{2}];

fig = figure;
scatter(X(:, 1), X(:, 2), '.')
axis image
min_white_margin(gca);
saveas(fig, 'data', 'png')

disp('Conducting some initial computations...')
% init comp
bound_x = [0 1];
bound_y = [0 max(X(:, 2))];
n = size(X, 1);
count = ones(n, 1);
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, bound_x, bound_y, count);
adj_mat = get_adj_mat( E, n );

fig = figure;
triplot(DT, 'Color', GRAY)
hold on
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
on_margin = find(cell_area> mean(cell_area(valid))+5*(var(cell_area(valid)).^(1/2)));
valid = setdiff(valid,on_margin);
invalid = [invalid; on_margin];

scatter(cx(valid), cy(valid), 12, cell_log_intensity(valid), 'filled')
colorbar
colormap(jet)
axis image
min_white_margin(gca);
saveas(fig, 'log_intensity', 'png')

disp('Getting initial seeds...')
% get seeds
st_x = 0.1;
en_x = 0.9;
st_y = 0.1*bound_y(2);
en_y = 0.9*bound_y(2);
step_x = 0.1;
step_y = 0.1*bound_y(2);
set_size = 20;
factor = 2;
k = 100;
set_size2 = 20;
[seeds, seeds_rej, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(st_x, en_x, st_y, en_y,...
   step_x, step_y, set_size, cell_log_intensity, cell_area, cx, cy, factor, k, set_size2, invalid,adj_mat);
num = num_s+num_s_pt

%  step_x = 0.8/63;
% step_y = 0.6/63;
%  [seeds_all, seeds_rej,num_s] = get_seeds_sim(st_x, en_x, st_y, en_y,...
%     step_x, step_y, set_size, cell_log_intensity, cell_area, cx, cy, factor,invalid,adj_mat);
% num = num_s;
 disp(['Number of regions is ', num2str(num)])

%plot the seeds
fig = figure;
triplot(DT, 'Color', GRAY)
hold on
% specify the colormap
colors = lines(num_s);
for i = 1:num_s
    scatter(cx(seeds{i}), cy(seeds{i}), 12, colors(i, :), 's', 'filled')
end
for i = 1:num_s_pt
    scatter(cx(seeds_pt{i}), cy(seeds_pt{i}), 12, 'r', 'd', 'filled')
end
for i = 1:length(seeds_rej)
    scatter(cx(seeds_rej{i}), cy(seeds_rej{i}), 12, 'k', '^', 'filled')
end
axis image
saveas(fig, 'seeds', 'png')

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

disp('Seeded region growing...')
% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid', true, 1000);
%disp(strcat("The time used for region growing is",num2str(toc)))
% plot the over-segmented image
% fig = figure;
% triplot(DT, 'Color', GRAY)
% hold on
% for i = 1:num_s
%     scatter(cx(region_sets{i}), cy(region_sets{i}), 12,  colors(i, :), 'filled')
% end
% for i = 1:num_s_pt
%     scatter(cx(region_sets{i+num_s}), cy(region_sets{i+num_s}), 12, 'r', 'filled')
% end
% axis image
% saveas(fig, 'over_seg', 'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('Region merging...')
[region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);
[sets_all, log_like_all] = merge_region_fast(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n);

%figure
%plot(log_like_all, '-o')

BIC_all = -2*log_like_all+6*(num-1:-1:0)'*log(n);
[min_BIC, index_BIC] = min(BIC_all);
elapsed_time = toc;
disp(strcat("The time used for region merging is",num2str(toc)))




fig = figure;
triplot(DT, 'Color', GRAY)
hold on
% the final result
for i = 1:num
    if ~isempty(selected{i})
        log_int = log(sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i}))/sum(cell_area(selected{i})));
        scatter(cx(selected{i}), cy(selected{i}), 12, log_int*ones(length(selected{i}), 1), 'filled')
    end
end
colorbar('SouthOutside')
colormap(hsv)
axis image
saveas(fig, 'final_seg', 'png')

disp(['Elapsed time is ', num2str(elapsed_time), ' seconds.'])

% save all variables (exclude figure handle)
clear fig
filename = 'real_full_fast_result';
current_datetime = clock;
for i = current_datetime
    filename = [filename, '_', num2str(fix(i))];
end
filename = [filename, '.mat'];
save(filename)
