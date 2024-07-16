

function [fourier_coef, contourx,contoury, area_poly,area_aic,area_bic, num_B_poly, num_B_aic,num_B_bic, intensity_B_poly,...
    intensity_B_aic, intensity_B_bic, num_O_poly,num_O_aic,num_O_bic, intensity_O_poly, intensity_O_aic,intensity_O_bic,...
    centx_poly, centy_poly, centx_aic, centy_aic,centx_bic, centy_bic] = sim_fit_random_real_data_withseed(X, seed, show_plot, P, threshold, rand_num, rep_itr, X_obs, min_select_AIC,min_select_BIC, penalty, M, seeds)
area_poly = 0; 
area_aic = 0;
num_B_poly = 0; 
num_B_aic = 0;
intensity_B_poly = 0;
intensity_B_aic = 0; 
num_O_poly = 0; 
num_O_aic = 0; 
intensity_O_poly = 0; 
intensity_O_aic = 0;
centx_poly = 0; 
centy_poly = 0; 
centx_aic = 0; 
centy_aic = 0;
fourier_coef = zeros(1,300);
contourx = [];
contoury = [];
% fit the model in simulation
rng(seed)
% init comp
[C,ix,~] = unique(X,'rows');
seeds = cell2mat(seeds);
seeds_indict = zeros(1,length(X));
seeds_indict(seeds) = 1;
X = C;
seeds_indict = seeds_indict(ix);

seeds = find(seeds_indict);
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

[invalid, ~] = get_invalid_cells(cell_log_intensity, adj_mat, n);
seeds = setdiff(seeds,invalid); %remove invalid cells from seeds
seeds = num2cell(seeds);

num = length(seeds);
disp(['Number of regions is ', num2str(num)])

% make a copy of variable seeds
region_sets = seeds;

% graph-based SRG
[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');

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
if n_region==1
    return 
else
    if n_region==2
        selected_post = selected_nonempty;
    else
    [region_intensity, region_area, region_num_cells, adj_mat_region] =...
        get_region_int_connect(n_region, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
    [~,~,sets_all_post, log_like_all_post,~,~,~] = merge_region_fast_beamsearch_optimization(n_region, region_area, ...
        region_intensity, selected_nonempty, adj_mat_region, region_num_cells, n, 20);
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
        %get aic result
        [invx,invy] = iFD(fourier_coef.',min_select_AIC);
        centx_aic = real(fourier_coef(1));
        centy_aic = imag(fourier_coef(1));
        area_aic = polyarea(invx,invy);
        num_B_aic = sum(inpolygon(X(:,1),X(:,2),invx,invy));
        intensity_B_aic  = num_B_aic/area_aic;
        num_O_aic = sum(inpolygon(X_obs(:,1),X_obs(:,2),invx,invy));
        intensity_O_aic = num_O_aic/area_aic;

        %get bic result
        [invx,invy] = iFD(fourier_coef.',min_select_BIC);
        centx_bic = real(fourier_coef(1));
        centy_bic = imag(fourier_coef(1));
        area_bic = polyarea(invx,invy);
        num_B_bic = sum(inpolygon(X(:,1),X(:,2),invx,invy));
        intensity_B_bic  = num_B_bic/area_bic;
        num_O_bic = sum(inpolygon(X_obs(:,1),X_obs(:,2),invx,invy));
        intensity_O_bic = num_O_bic/area_bic;
        
    end
end

if show_plot
         GRAY = [0.6 0.6 0.6];

%      figure
%      colors = lines(n_region);
%     triplot(DT, 'Color', GRAY)
%     hold on
%     index = 0;
%     for i = 1:n_region
%         if ~isempty(selected_nonempty{i})
%             index = index+1;
%             scatter(cx(selected_nonempty{i}), cy(selected_nonempty{i}), 12, colors(index, :), 'filled')
%         end
%     end
%     axis image
%     figure
%     GRAY = [0.6 0.6 0.6];
%     colors = lines(n_region);
%     triplot(DT, 'Color', GRAY)
%     hold on
%     index = 0;
%     for i = 1:n_region
%         if ~isempty(selected_post{i})
%             index = index+1;
%             scatter(cx(selected_post{i}), cy(selected_post{i}), 12, colors(index, :), 'filled')
%         end
%     end
%     axis image
%     figure 
    hold on
    scatter(X_obs(:,1),X_obs(:,2),2,GRAY,'filled')
    axis equal
    axis([0,1,0,1])
    plot(x_sample,y_sample,'Color','r','LineWidth',2)

end

end
