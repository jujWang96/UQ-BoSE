close all
clear
GRAY = [0.6,0.6,0.6];

cd('~/Desktop/gsrg/Astro_sim/real_img/Hamburger')
addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))
fid=fopen('Hamburger_11007/hamburger_0.txt');
C=textscan(fid,'%f%f%f%f%f%f%f%f%f','HeaderLines',11,'delimiter','\r\n','CommentStyle','*');
fclose(fid);
filename = 'real_0_grid_result';
penalty_term = 6;
filename = strcat(filename,'_penalty',num2str(penalty_term));
int_factor = 2; %the intensity factor to scale the partial segmentation to full segmentation 

X = [C{1},C{2}];

fig = figure;
scatter(X(:, 1), X(:, 2), '.')
axis image
min_white_margin(gca);
saveas(fig, 'data', 'epsc')

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
step_x = 0.2;
step_y = 0.2*bound_y(2);
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
colors = lines(num);
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
%plot the over-segmented image
fig = figure;
plot(X(valid,1),X(valid,2),'.')
hold on
growed = 0;
for i = 1:num_s
    scatter(cx(region_sets{i}), cy(region_sets{i}), 12,  colors(i, :))
    growed = growed+ length(region_sets{i});
end
for i = 1:num_s_pt
    scatter(cx(region_sets{i+num_s}), cy(region_sets{i+num_s}), 12, 'r')
    growed = growed+ length(region_sets{i+num_s});

end
growed
axis image
% saveas(fig, 'over_seg', 'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rep_itr = 30000;
rand_num = 3;
disp('Region merging...')
[region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);
% [sets_all, log_like_all] = merge_region_fast(num, region_area, ...
%     region_intensity, region_sets, adj_mat_region, region_num_cells, n);


selected_rep = cell(1,rep_itr);
BIC_rep = zeros(1,rep_itr);
tic;
%repeat multiple times
parfor ii = 1:rep_itr
    rng(ii)
    [sets_all, log_like_all] = merge_region_random_fast(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n, rand_num);
    BIC_all = -2*log_like_all+penalty_term*(num-1:-1:0)'*log(n);
    [min_BIC, index_BIC] = min(BIC_all);
    selected = sets_all{index_BIC};
    selected_rep{ii} = selected;
    BIC_rep(ii) = min_BIC;
end

[min_BIC, index_rep] = min(BIC_rep);
selected = selected_rep{index_rep};
elapsed_time = toc;



fig = figure;

triplot(DT, 'Color', GRAY)
hold on
% the final result
for i = 1:num
    if ~isempty(selected{i})
        log_int = log(sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i})*int_factor)/sum(cell_area(selected{i})));
        scatter(cx(selected{i}), cy(selected{i}), 12, log_int*ones(length(selected{i}), 1), 'filled')
    end
end
colorbar('SouthOutside')
colormap(hsv)
axis image
saveas(fig, 'final_seg', 'epsc')

disp(strcat("The minimum BIC value is ", num2str(min_BIC)))
figure;
plot_segmentation_wo_voronoi(DT, selected, cx, cy,colors ,8,true)

disp(['Elapsed time is ', num2str(elapsed_time), ' seconds.'])

% save all variables (exclude figure handle)
clear fig
current_datetime = clock;
for i = current_datetime
    filename = [filename, '_', num2str(fix(i))];
end
filename_mat = [filename, '.mat'];
save(filename_mat)

%% write to txt file 
index = 0;
photon_id = [];
region_id = [];
for i = 1:length(selected)
    if ~isempty(selected{i})
        index = index + 1;
        photon_id = [photon_id,selected{i}];
        region_id = [region_id,ones(1,length(selected{i}))*index];
    end
end
%include invalid
photon_id = [photon_id,invalid'];
region_id = [region_id,ones(1, length(invalid))*(-1)];
[~,idx]=sort(photon_id);

region_id = region_id(idx);
fildId = fopen([filename,'.txt'],'w');
formatSpec = '%f    %f    %d\n';
fprintf(fildId,formatSpec,[X(:,1)';X(:,2)';region_id]);


