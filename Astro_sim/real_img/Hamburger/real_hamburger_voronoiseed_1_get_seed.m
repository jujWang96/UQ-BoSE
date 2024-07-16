close all
clear
GRAY = [0.6,0.6,0.6];
cd('~/Desktop/gsrg/Astro_sim/real_img/Hamburger')
addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))
addpath(genpath('~/Desktop/gsrg/plot_util'))

fid=fopen('Hamburger_11007/hamburger_1.txt');
C=textscan(fid,'%f%f%f%f%f%f%f%f%f','HeaderLines',11,'delimiter','\r\n','CommentStyle','*');
fclose(fid);
filename = 'real_1_voronoi_result';
penalty_term = 6;
filename = strcat(filename,'_penalty',num2str(penalty_term));
M = 150; %run greedy merge until 150 region left
X = [C{1},C{2}];

fig = figure;
scatter(X(:, 1), X(:, 2), '.')
axis image
min_white_margin(gca);
saveas(fig, 'data', 'epsc')

disp('Conducting some initial computations...')
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
P=15;
threshold = 8;
 [seeds, num_region] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat,cell_area,P, threshold);
num = num_region;
%  step_x = 0.8/63;
% step_y = 0.6/63;
%  [seeds_all, seeds_rej,num_s] = get_seeds_sim(st_x, en_x, st_y, en_y,...
%     step_x, step_y, set_size, cell_log_intensity, cell_area, cx, cy, factor,invalid,adj_mat);
% num = num_s;
 disp(['Number of regions is ', num2str(num)])

%plot the seeds
fig = figure;
colors = lines(num);
plot_seeds_one(DT, cx, cy, seeds, [], [], colors, num, 0)
set(gca, 'fontsize', 18)
min_white_margin(gca);
set(fig, 'Position', [0, 0, 490, 380]);

saveas(fig, 'seeds', 'png')

