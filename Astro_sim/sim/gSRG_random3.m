function [seg_contour,seg_region_intensity,seg_selected,seg_min_BIC] = gSRG_random3(X,draw_s,draw_c, Color,N_r,bound,drop,rand_num,rand_itr,itr_num,m) 

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
%
%Output:
% C: cell array, each element contains x1,x2,y1,y2 representing the 
% contour of of one region 

subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);


% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
[seeds, seeds_rej, seeds_pt, num_s, num_s_pt] = get_seeds_sim_local_max(0.1, 0.9, 0.1, 0.9,...
    0.1, 0.1, 5, cell_log_intensity, cell_area, cx, cy, 2, 30, 5, invalid,adj_mat);
num = num_s+num_s_pt;
frame = [0.45 0.55 0.45 0.55];
pt_size = 0.2;
% plot the seeds

seeds_all = [seeds seeds_pt];

% make a copy of variable seeds
region_sets = seeds_all;

% graph-based SRG

[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');
glb_min_BIC = realmax;
[sets_all_def, log_like_all_def,region_intensity_def,region_area,region_num_cells,region_indices,adj_mat_region] = merge_def(num, cell_area, ...
        cell_log_intensity, region_sets, adj_mat, n,rand_itr);
glb_min_BIC = realmax;
seg_min_BIC = realmax;
for itr = 1:itr_num   
    [sets_all, log_like_all,region_intensity_all] = merge_rand(num, sets_all_def, log_like_all_def,region_intensity_def,region_area,region_num_cells,cell_area, ...
    cell_log_intensity, region_sets, adj_mat_region, n,region_indices,rand_num,rand_itr);
    BIC_all = -2*log_like_all+m*(num-1:-1:0)'*log(n);
    [min_BIC, min_index] = min(BIC_all);
    
    %min_index = get_min_index(N_r,bound,BIC_all);
    %min_BIC = BIC_all(min_index);
   % [min_BIC,min_index] = min(BIC_all);
    if min_BIC<glb_min_BIC
        glb_min_BIC = min_BIC;
        glb_selected = sets_all{min_index};
        glb_region_intensity = region_intensity_all{min_index};
    end
    if BIC_all(length(BIC_all)-1)<seg_min_BIC
        seg_min_BIC = BIC_all(length(BIC_all)-1);
        seg_selected = sets_all{length(BIC_all)-1};
        seg_region_intensity = region_intensity_all{length(BIC_all)-1};
   end
end
% % 
[glb_contour,~] = get_contour(n, DT,glb_selected,adj_mat,invalid);
[seg_contour,~] = get_contour(n, DT,seg_selected,adj_mat,invalid);

%subplot(2,2,1)
%plot_segmentation_pause(DT, glb_selected, cx, cy, lines(num),0,glb_contour,frame,0.3)
%subplot(2,2,2)
%plot_intensity( cx,cy, glb_selected, glb_contour,glb_region_intensity,0.3)
%axis(frame)
%leg = legend(string(glb_contour.regionId))
%leg.FontSize = 5;
%set(gcf, 'Position', [50 50 600 600]); %

%subplot(2,2,3)
%plot_segmentation_pause(DT, seg_selected, cx, cy, lines(num),0,seg_contour,frame,0.3)

%function min_index = get_min_index(N_r,bound,BIC_all)
%	L = length(BIC_all); 
%	if isempty(N_r)
 %               disp('no bound')
  %  		[~, min_index] = min(BIC_all); 
%	else
 %   		if ~bound && N_r<=L
  %      		min_index = L+1-N_r;
   % 		elseif bound && N_r<=L
    %    		[~, min_index] = min(BIC_all(L+1-N_r:end));
     %   		min_index = min_index+L-N_r;
    %		end
%	end

%end

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

    function [sets_all_def, log_like_all_def,region_intensity_all] = merge_rand(num, sets_all_def, log_like_all_def,region_intensity,region_area,region_num_cells,cell_area,...
    cell_log_intensity, region_sets, adj_mat_region, n,region_indices,rand_num,rand_itr)
        for t = num-rand_itr:num-1
            max_incs = -realmax*ones(1,rand_num);
            merge_regions1 = zeros(1,rand_num);
            merge_regions2 = zeros(1,rand_num);
            reg_num = 0;
            for i = region_indices
                for j = region_indices
                    if i<j && adj_mat_region(i, j)
                        %disp('there is region connected')
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
            idx = randi(min(rand_num,reg_num));

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


