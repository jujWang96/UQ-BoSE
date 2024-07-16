close all
clear
GRAY = [0.6,0.6,0.6];
cd('~/Desktop/gsrg/Astro_sim/real_img/G292')
addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))
addpath(genpath('~/Desktop/gsrg/plot_util'))
filename = 'hard_12555';
fid=fopen(strcat('SNR_G292.0+1.8_photons/',filename,'.txt'));
C=textscan(fid,'%f%f%f%f%f%f%f%f','HeaderLines',12,'delimiter','\r\n','CommentStyle','*');
fclose(fid);
%filename = 'real_full_voronoi_result';
filename = strcat(filename, '_result_');

penalty_term = 6;
%filename = strcat(filename,'_penalty',num2str(penalty_term));
M = 150; %run greedy merge until 150 region left
X = [C{1},C{2}];

fig = figure;
scatter(X(:, 1), X(:, 2), '.')
axis image
min_white_margin(gca);
%saveas(fig, 'data', 'epsc')

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
figure
scatter(cx(valid), cy(valid), 12, cell_log_intensity(valid), 'filled')
colorbar
colormap(jet)
axis image
hold on
scatter(cx(on_margin),cy(on_margin),20, 'k')
min_white_margin(gca);
%saveas(fig, 'log_intensity', 'png')

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
P=20;
threshold = 5;
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
triplot(DT, 'Color', GRAY)
hold on
% specify the colormap
colors = lines(num);
for i = 1:num
    scatter(cx(seeds{i}), cy(seeds{i}), 12, colors(i, :), 's', 'filled')
end

axis image
%saveas(fig, 'seeds', 'png')

seeds_all = seeds;

% make a copy of variable seeds
region_sets = seeds_all;

disp('Seeded region growing...')
% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid', true, 1000);
figure;
plot_segmentation_wo_voronoi(DT, region_sets, cx, cy,colors ,8,true)

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
%greedy merge until M region left then do randomized search
if num> M
    [sets_all_greedy, log_like_all_greedy]  = merge_region_fast(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n);
    sets_greedy = sets_all_greedy{num-M+1};
    region_sets = sets_greedy(~cellfun(@isempty,sets_greedy));
    num = length(region_sets); %should be equal to M
    [region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);
end

figure;
plot_segmentation_wo_voronoi(DT, region_sets, cx, cy,colors ,8,true)


%plot sets_all_greedy last 
sets_debug = sets_all_greedy{num};
plot_segmentation_wo_voronoi(DT, sets_debug, cx, cy,colors ,8,true)

rep_itr = 10000;
rand_num = 5;
selected_rep = cell(1,rep_itr);
BIC_rep = zeros(1,rep_itr);
tic;
%repeat multiple times


%debug
% [sets_all_greedy, log_like_all_greedy]  = merge_region_fast(num, region_area, ...
%     region_intensity, region_sets, adj_mat_region, region_num_cells, n);
% sets_debug = sets_all_greedy{num};
% figure
% plot_segmentation_wo_voronoi(DT, sets_debug, cx, cy,colors ,8,true)
% 
% [sets_all, log_like_all] = merge_region_random_fast_debug(num, region_area, ...
%     region_intensity, region_sets, adj_mat_region, region_num_cells, n, 1,ii, DT,cx,cy);

%%

parfor ii = 1:rep_itr
    [sets_all, log_like_all] = merge_region_random_fast(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n, rand_num,ii);
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
        log_int = log(sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i}))/sum(cell_area(selected{i})));
          scatter(cx(selected{i}), cy(selected{i}), 12, log_int*ones(length(selected{i}), 1), 'filled')
    end
end
colorbar('SouthOutside')
colormap(hsv)
axis image
%saveas(fig, 'final_seg', 'epsc')
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
segment_intensity = [];
segment_area = [];
for i = 1:length(selected)
    if ~isempty(selected{i})
        index = index + 1;
        photon_id = [photon_id,selected{i}];
        region_id = [region_id,ones(1,length(selected{i}))*index];
        segment_area = [segment_area,ones(1,length(selected{i}))*sum(cell_area(selected{i}))];
        segment_intensity = [segment_intensity,ones(1,length(selected{i}))*sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i}))/sum(cell_area(selected{i}))];
      
    end
end
%include the isolated cells to invalid
isolate_id = setdiff(valid, photon_id);
invalid = [invalid;isolate_id];
%include the invalid cells
photon_id = [photon_id,invalid'];
region_id = [region_id,ones(1, length(invalid))*(-1)];
segment_area = [segment_area, zeros(1, length(invalid))];
segment_intensity = [segment_intensity, zeros(1, length(invalid))];

[~,idx]=sort(photon_id);

region_id = region_id(idx);
segment_area = segment_area(idx);
segment_intensity = segment_intensity(idx);
fildId = fopen([filename,'.txt'],'w');
formatSpec = '%f    %f   %f    %f    %f    %f     %f    %f    %d   %f    %f    %f \n';
fprintf(fildId,formatSpec,[C{1}';C{2}';C{3}';C{4}';C{5}';C{6}'; ...
    C{7}';C{8}';region_id;segment_intensity;segment_area;cell_area']);



%debug 
figure
scatter(C{1},C{2},12,region_id);