

function [fourier_coef, contourx,contoury, area_poly,area_fd, num_B_poly, num_B_fd, intensity_B_poly,...
    intensity_B_fd, num_O_poly,num_O_fd, intensity_O_poly, intensity_O_fd,...
    centx_poly, centy_poly, centx_fd, centy_fd] = sim_fit_2300chandra_bootstrap(X, seed, show_plot, P, threshold, rand_num, rep_itr, X_obs, min_modelselect,penalty, M,origin_invx,origin_invy,bootmethod)
X = unique(X,'rows');
X_obs = unique(X_obs,'rows');
area_poly = 0; 
area_fd = 0;
num_B_poly = 0; 
num_B_fd = 0;
intensity_B_poly = 0;
intensity_B_fd = 0; 
num_O_poly = 0; 
num_O_fd = 0; 
intensity_O_poly = 0; 
intensity_O_fd = 0;
centx_poly = 0; 
centy_poly = 0; 
centx_fd = 0; 
centy_fd = 0;
fourier_coef = zeros(1,300);
contourx = [];
contoury = [];
% fit the model in simulation
rng(seed)
imagename = '2300chandra_vor';
imagename = strcat(imagename,num2str(bootmethod),'_');
imagename = strcat(imagename, 'P-',num2str(P),'S-',num2str(threshold),'p-',num2str(penalty) ,'seed-',num2str(seed));
% init comp
figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);
ncol = 4;
n_plot = 10;
idx_plot = 0;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
[~, ~, ~, ~, ~, ~, cell_area_obs] = init_comp(X_obs, [0 1], [0 1], ones(size(X_obs, 1), 1));

adj_mat = get_adj_mat( E, n );
idx_plot = idx_plot+1;
subplot(ceil(n_plot/ncol),ncol,idx_plot)

histogram(cell_area_obs, "FaceAlpha", 0.5)
hold on
histogram(cell_area,"FaceAlpha", 0.5)

% get seeds
[invalid, ~] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[seeds, num_region] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat,cell_area,P, threshold);
num = num_region;
%disp(['Number of regions is ', num2str(num)])


idx_plot = idx_plot+1;

gca = subplot(ceil(n_plot/ncol),ncol,idx_plot);
plot_seeds(DT, cx, cy, seeds, [], [], lines(num), num, 0);
gca = image_rescaling(gca,'NGC2300');

% make a copy of variable seeds
region_sets = seeds;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');
idx_plot = idx_plot+1;
subplot(ceil(n_plot/ncol),ncol,idx_plot)
plot_segmentation_wo_voronoi(DT, region_sets, cx, cy,lines(num) ,8,true);
gca = image_rescaling(gca,'NGC2300');
% get source boundaries formed by Voronoi cell edges
hold on
[V, R] = voronoiDiagram(DT);

vx_edges_all = {};
vy_edges_all = {};
for j = 1:num
    cells_in_region = R(region_sets{j});
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
for j = 1:num
    % original red is brighter
    line(vx_edges_all{j}, vy_edges_all{j}, 'Color', 'k', 'LineWidth',3)
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

if num> M
    [sets_all_greedy, log_like_all_greedy]  = merge_region_fast(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n);
    sets_greedy = sets_all_greedy{num-M+1};
    region_sets = sets_greedy(~cellfun(@isempty,sets_greedy));
    num = length(region_sets); %should be equal to M
    [region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);
end
selected_rep = cell(1,rep_itr);
BIC_rep = zeros(1,rep_itr);

%repeat multiple times
parfor ii = 1:rep_itr
    [sets_all, log_like_all] = merge_region_random_fast(num, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n,rand_num, ii);
    BIC_all = -2*log_like_all+penalty*(num-1:-1:0)'*log(n);
    [min_BIC, index_BIC] = min(BIC_all);
    selected = sets_all{index_BIC};
    selected_rep{ii} = selected;
    BIC_rep(ii) = min_BIC;
end

[min_BIC, index_rep] = min(BIC_rep);
selected = selected_rep{index_rep};
idx_plot = idx_plot+1;
subplot(ceil(n_plot/ncol),ncol,idx_plot)
plot_segmentation_wo_voronoi(DT, selected, cx, cy,lines(num) ,8,true);
gca = image_rescaling(gca,'NGC2300');


n_region = 0;
pred_class_all = zeros(n, 1);
for i = 1:num
    if ~isempty(selected{i})
        n_region = n_region+1;
        pred_class_all(selected{i}) = n_region;
    end
end


%do another beam merge until two regions left
selected_nonempty = selected(~cellfun(@isempty,selected));
idx_plot = idx_plot+1;
subplot(ceil(n_plot/ncol),ncol,idx_plot)
%plot the area 

hold on
colors = lines(n_region);
area_regions = [];
for i =1:n_region
    barh(i, 1, FaceColor =colors(i,:))
    area_regions = [area_regions,sum(cell_area(selected_nonempty{i}))];
end
convertrate = 34.7*33.5*100;
fs = 15;
text(ones(1,n_region)*1.5,(1:n_region),num2str(convertrate*area_regions'),'vert','middle','horiz','center','FontSize',fs); 
box off
xlim([0,2])

idx_plot = idx_plot+1;
subplot(ceil(n_plot/ncol),ncol,idx_plot)
hold on
colors = lines(n_region);
flux_regions = [];
for i =1:n_region
    barh(i, 1, FaceColor =colors(i,:))
    flux_regions = [flux_regions,length((selected_nonempty{i}))];
end
text(ones(1,n_region)*1.5,(1:n_region),num2str(flux_regions'),'vert','middle','horiz','center','FontSize',fs); 
box off
xlim([0,2])

if n_region==1
    return 
else
    if n_region==2
        selected_post = selected_nonempty;
    else
    [region_intensity, region_area, region_num_cells, adj_mat_region] =...
        get_region_int_connect(n_region, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
    [~,~,sets_all_post, log_like_all_post,~,~,~] = merge_region_fast_beamsearch_optimization(n_region, region_area, ...
        region_intensity, selected_nonempty, adj_mat_region, region_num_cells, n, 200);
    selected_post = sets_all_post{n_region-1};
    end
    %get the countour and calculate fd
    [seg_contour,~] = get_contour(n, DT,selected_post,adj_mat,invalid,cell_area);
    objId = setdiff(seg_contour.regionId,seg_contour.background);

    if length(objId)>1 || length(objId)<1

        return;
    else
        [contourx,contoury] = get_curve(seg_contour.contourV{objId},false);
        area_poly = polyarea(contourx,contoury);
        num_B_poly = sum(inpolygon(X(:,1),X(:,2),contourx,contoury));
        intensity_B_poly  = num_B_poly/area_poly;
        num_O_poly = sum(inpolygon(X_obs(:,1),X_obs(:,2),contourx,contoury));
        intensity_O_poly = num_O_poly/area_poly;
        %smooth curve 
        [centx_poly,centy_poly] = centroid(polyshape(contourx,contoury));
        [x_sample,y_sample] = sample_curve(contourx,contoury,300,centx_poly, false);
        fourier_coef = FD(x_sample,y_sample).';
        [invx,invy] = iFD(fourier_coef.',min_modelselect);
        centx_fd = real(fourier_coef(1));
        centy_fd = imag(fourier_coef(1));
        area_fd = polyarea(invx,invy);
        num_B_fd = sum(inpolygon(X(:,1),X(:,2),invx,invy));
        intensity_B_fd  = num_B_fd/area_fd;
        num_O_fd = sum(inpolygon(X_obs(:,1),X_obs(:,2),invx,invy));
        intensity_O_fd = num_O_fd/area_fd;
    end
end
idx_plot = idx_plot+1;
subplot(ceil(n_plot/ncol),ncol,idx_plot)
plot_segmentation_wo_voronoi(DT, selected_post, cx, cy,lines(num) ,8,true);
gca = image_rescaling(gca,'NGC2300');





idx_plot = idx_plot+1;
subplot(ceil(n_plot/ncol),ncol,idx_plot)

gca = scatter(X(:,1),X(:,2),1,'k.');
hold on 
plot(x_sample,y_sample,'Color','r','LineWidth',2)
plot(origin_invx,origin_invy,'Color','b','LineWidth',2)
gca = image_rescaling(gca,'NGC2300');

axis square


selected_nonempty_post = selected_post(~cellfun(@isempty,selected_post));

idx_plot = idx_plot+1;
subplot(ceil(n_plot/ncol),ncol,idx_plot)
hold on
colors = lines(2);
area_regions = [];
for i =1:2
    barh(i, 1, FaceColor =colors(i,:))
    area_regions = [area_regions,sum(cell_area(selected_nonempty_post{i}))];
end
convertrate = 34.7*33.5*100;

text(ones(1,2)*1.5,(1:2),num2str(convertrate*area_regions'),'vert','middle','horiz','center','FontSize',fs); 
box off
xlim([0,2])

idx_plot = idx_plot+1;
subplot(ceil(n_plot/ncol),ncol,idx_plot)
hold on
colors = lines(2);
flux_regions = [];
for i =1:2
    barh(i, 1, FaceColor =colors(i,:))
    flux_regions = [flux_regions,length((selected_nonempty_post{i}))];
end
text(ones(1,2)*1.5,(1:2),num2str(flux_regions'),'vert','middle','horiz','center','FontSize',fs); 
box off
xlim([0,2])



set(gcf, 'Position', [50 50 1300 800]); %

saveas(gcf, strcat(imagename),'png')

end
