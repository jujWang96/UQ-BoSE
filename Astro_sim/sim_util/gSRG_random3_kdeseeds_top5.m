function [contours,glb_region_intensitys,selecteds,min_BICs,seg_nums,glb_min_BICs] = gSRG_random3_kdeseeds_top5(X,draw_s,draw_c, Color,N_r,bound,drop,rand_num,rand_itr,itr_num,m,threshold) 

subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.02], [0.02 0.02], [0.02 0.02]);

%
%
%Use SRG-graph to segment the points, merge them uniform randomly and
%return the contour of the segmentation. 
%Input: 
% X: all of the points
% draw_s: logic; if true then draw the segmentation
% draw_c: logic; if true then draw the contour
% Color: the color of contour
% N_r: the true number of region 
% drop: logic, if true then drop this simulation when the segmentation 
% not equal to true numbe of region
% restrict: if the segmentation region number is less than restrict number, 
% return segmentation of N_r regions 
% rand_num: number of possible merge to choose each time
% rand_itr: the iteration that uses random merge
% m: penalty coefficient
% seeds_all: self-defined seeds.
%Output:
% C: cell array, each element contains x1,x2,y1,y2 representing the 
% contour of of one region 



% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
%[seeds, num] = get_seeds_sim_kde(cx, cy,invalid, adj_mat, threshold);
[seeds_all, num] = get_seeds_sim_kde(cx, cy,invalid, adj_mat, 5);
num
% plot the seeds
colors = lines(num);
region_sets = seeds_all;

% graph-based SRG

[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');
[sets_all_def, log_like_all_def,region_intensity_def,region_area,region_num_cells,region_indices,adj_mat_region] = merge_def(num, cell_area, ...
        cell_log_intensity, region_sets, adj_mat, n,rand_itr);
glb_min_BICs = realmax*ones(1,5);
min_BICs = zeros(1,itr_num);
seg_nums = zeros(1,itr_num);
for itr = 1:itr_num   
    [sets_all, log_like_all,region_intensity_all] = merge_rand(num, sets_all_def, log_like_all_def,region_intensity_def,region_area,region_num_cells,cell_area, ...
    cell_log_intensity, region_sets, adj_mat_region, n,region_indices,rand_num,rand_itr);
    BIC_all = -2*log_like_all+m*(num-1:-1:0)'*log(n);
    [min_BIC, min_index] = min(BIC_all);
[max_glb_BIC,max_glb_index] = max(glb_min_BICs);
    if min_BIC<max_glb_BIC && ~any([glb_min_BICs]==min_BIC)
        glb_min_BICs(max_glb_index) = min_BIC;
        selecteds{max_glb_index} = sets_all{min_index};
        glb_region_intensitys{max_glb_index} = region_intensity_all{min_index};
    end
    min_BICs(itr) = min_BIC;
    seg_nums(itr) = length(BIC_all)-min_index+1;

    %min_index = get_min_index(N_r,bound,BIC_all);
    %min_BIC = BIC_all(min_index);
   % [min_BIC,min_index] = min(BIC_all);
end
%subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);
order = 0;
 [sort_glb_BIC,sort_glb_index] = sort(glb_min_BICs);

for output = sort_glb_index
selected = selecteds{output};
glb_region_intensity = glb_region_intensitys{output};
order = order+1;
pt_size = 0.05
subplot(2,2,1)
plot_seeds(DT, cx, cy, seeds_all, [], [], colors, num, 0)

subplot(2,2,2)
if draw_s ==true
    colors = lines(num);
    %colors = distinguishable_colors(num);
    plot_segmentation(DT, selected, cx, cy, colors,pt_size);
    axis([0 1 0 1])   
    axis square;
    hold on
end
[contour,N_s] = get_contour(n, DT, selected,adj_mat,invalid);
contours{output} = contour;
if draw_c == true  
        for seg = contour.regionId
            contourVcurr = contour.contourV{seg};
                x1 = contourVcurr(:,1);
                x2 = contourVcurr(:,2);
                y1 = contourVcurr(:,3);
                y2 = contourVcurr(:,4);
                plot([x1';x2'], [y1';y2'],'Color',Color,'LineWidth',0.1);
                axis([0 1 0 1])   
                axis square;
                hold on;
           end

end
subplot(2,2,4)
plot_intensity(cx,cy, selected, contour, glb_region_intensity,pt_size)
set(gcf, 'Position', [50 50 600 600]); %
saveas(gcf,strcat('Arp299_kdeseedsgt5_top',num2str(order),'_150.pdf'))
close
disp('plot updated')
end
%disp('plot intensity')
%set(gcf, 'Position', [0, 0, 1000, 560]);
%drop if the segmentation touches the boundary



function min_index = get_min_index(N_r,bound,BIC_all)
	L = length(BIC_all); 
	if isempty(N_r)
                disp('no bound')
    		[~, min_index] = min(BIC_all); 
	else
    		if ~bound && N_r<=L
        		min_index = L+1-N_r;
    		elseif bound && N_r<=L
        		[~, min_index] = min(BIC_all(L+1-N_r:end));
        		min_index = min_index+L-N_r;
    		end
	end

end

 function [sets_all_def, log_like_all_def,region_intensity,region_area,region_num_cells,region_indices,adj_mat_region] = merge_def(num, cell_area, ...
        cell_log_intensity, region_sets, adj_mat, n,rand_itr)
    [region_intensity, region_area, region_num_cells, adj_mat_region] =...
        get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);

    % compute the log-likelihood of the over-segmented graph before region merging
        log_like = get_metric_value('log-like', n, region_sets, region_area, region_intensity, region_num_cells);
        % construct a region merging tree
        % t represents the iter # of region merging
        region_indices = 1:num;
        log_like_all_def = zeros(num, 1);
        log_like_all_def(1) = log_like;
        sets_all_def = cell(num, 1);
        sets_all_def{1} = region_sets;
        for t = 1:num-rand_itr-1
            % find the pair of regions with the maximal log-likelihood increasement
             max_inc= -realmax;
             for i = region_indices
                for j = region_indices
                    if i<j && adj_mat_region(i, j)
                        numerator = region_intensity(i)*region_area(i)+region_intensity(j)*region_area(j);
                        sum_area = region_area(i)+region_area(j);
                        delta_log_like = region_num_cells(i)*log(numerator/region_intensity(i)/sum_area)+...
                            region_num_cells(j)*log(numerator/region_intensity(j)/sum_area);
                        % delta_log_like is defined as new log-like minus old
                        % log-like
                        % update the max delta value 
                        if delta_log_like>max_inc
                            max_inc = delta_log_like;
                            merge_region1 = i;
                            merge_region2 = j;
                        end
                    end
                end
             end
            % del merge_region2 from region_indices since it will be merged with merge_region1
            region_indices(region_indices==merge_region2) = [];
            % update log-likelihood
            log_like_all_def(t+1) = log_like_all_def(t)+max_inc;
            % merge merge_region1 and merge_region2
            % only keep merge_region1
            region_intensity(merge_region1) = (region_intensity(merge_region1)*region_area(merge_region1)+...
                region_intensity(merge_region2)*region_area(merge_region2))/...
                (region_area(merge_region1)+region_area(merge_region2));
            region_area(merge_region1) = region_area(merge_region1)+region_area(merge_region2);
            region_num_cells(merge_region1) = region_num_cells(merge_region1)+region_num_cells(merge_region2);
            % update the connectivity
            idx = setdiff(find(adj_mat_region(merge_region2, :)==1), merge_region1);
            adj_mat_region(merge_region1, idx) = 1;
            adj_mat_region(idx, merge_region1) = 1;
            % record the sets after region merging
            tmp = sets_all_def{t};
            tmp{merge_region1} = [tmp{merge_region1} tmp{merge_region2}];
            tmp{merge_region2} = [];
            sets_all_def{t+1} = tmp;
        end


    end

    function [sets_all_def, log_like_all_def,region_intensity_all] = merge_rand(num, sets_all_def, log_like_all_def,region_intensity,region_area,region_num_cells,cell_area, ...
    cell_log_intensity, region_sets, adj_mat_region, n,region_indices,rand_num,rand_itr)
        for t = num-rand_itr:num-1
            max_incs = -realmax*ones(1,rand_num);
            merge_regions1 = zeros(1,rand_num);
            merge_regions2 = zeros(1,rand_num);
            reg_num = 0;
            for i = region_indices
                for j = region_indices
                    if i<j && adj_mat_region(i, j)
                        numerator = region_intensity(i)*region_area(i)+region_intensity(j)*region_area(j);
                        sum_area = region_area(i)+region_area(j);
                        delta_log_like = region_num_cells(i)*log(numerator/region_intensity(i)/sum_area)+...
                            region_num_cells(j)*log(numerator/region_intensity(j)/sum_area);
                   % delta_log_like is defined as new log-like minus old
                        % log-like
                        % update the max delta value 
                        [M,I] = min(max_incs);
                        if delta_log_like>M
                           reg_num = reg_num+1;
                            max_incs(I) = delta_log_like;
                            merge_regions1(I) = i;
                            merge_regions2(I) = j;
                        end
                    end
                end
            end
            %randomly sample the merge based on probability distribution of softmax transformation value of delta_log_like
            idx = randi(min(rand_num,reg_num));
            %idx = randi(1:min(rand_num,reg_num))

            merge_region1=merge_regions1(idx);
            merge_region2= merge_regions2(idx);
            max_inc = max_incs(idx);
            % del merge_region2 from region_indices since it will be merged with merge_region1
            region_indices(region_indices==merge_region2) = [];
            % update log-likelihood
            log_like_all_def(t+1) = log_like_all_def(t)+max_inc;
            % merge merge_region1 and merge_region2
            % only keep merge_region1
            region_intensity(merge_region1) = (region_intensity(merge_region1)*region_area(merge_region1)+...
                region_intensity(merge_region2)*region_area(merge_region2))/...
                (region_area(merge_region1)+region_area(merge_region2));
            region_area(merge_region1) = region_area(merge_region1)+region_area(merge_region2);
            region_num_cells(merge_region1) = region_num_cells(merge_region1)+region_num_cells(merge_region2);
            % update the connectivity
            index = setdiff(find(adj_mat_region(merge_region2, :)==1), merge_region1);
            adj_mat_region(merge_region1, index) = 1;
            adj_mat_region(index, merge_region1) = 1;
            % record the sets after region merging
            tmp = sets_all_def{t};
            tmp{merge_region1} = [tmp{merge_region1} tmp{merge_region2}];
            tmp{merge_region2} = [];
            sets_all_def{t+1} = tmp;
	    region_intensity_all{t+1} = region_intensity;

        end
        
        
    end
end


