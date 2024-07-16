function [contour,selected] = gSRG_random(X,draw_s,draw_c, Color,N_r,drop,restrict,rand_num,rand_itr,itr_num,m) 
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


% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );

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

[region_sets, labeled_cells] = SRG_graph(region_sets, cell_log_intensity, cell_area, n, adj_mat, invalid');
glb_min_BIC = realmax;
[sets_all_def, log_like_all_def,region_intensity,region_area,region_num_cells,region_indices,adj_mat_region] = merge_def(num, cell_area, ...
        cell_log_intensity, region_sets, adj_mat, n,rand_itr);
glb_min_BIC = realmax;
for itr = 1:itr_num   
    [sets_all, log_like_all] = merge_rand(num, sets_all_def, log_like_all_def,region_intensity,region_area,region_num_cells,cell_area, ...
    cell_log_intensity, region_sets, adj_mat_region, n,region_indices,rand_num,rand_itr);
    

BIC_all = -2*log_like_all+m*(num-1:-1:0)'*log(n);
    [min_BIC, min_index] = min(BIC_all);

    [L,~] = size(sets_all);
    N_s = L+1-min_index;
    %force to merge 
    if N_s <= restrict
        min_index = L+1-N_r;
        min_BIC = BIC_all(min_index);
    end

    if min_BIC<glb_min_BIC
        glb_min_BIC = min_BIC;
        selected = sets_all{min_index};
    end
end
disp(glb_min_BIC)

if draw_s == true
    figure
    colors = lines(num);
    plot_segmentation(DT, selected, cx, cy, colors);
    axis([0 1 0 1])   
    axis square;
disp('plotted')
end
%outputID = fopen(strcat('n132doutputg50000_BIC12/','s005_10000_gSRGoutput_r1_part',num2str(p),'.txt'),'w');
 %   for n = 1:length(X)
  %      for s =1:length(selected)
   %         seg = selected{s};
    %        if any(seg(:) == n)
     %           fprintf(outputID,'%f  %f  %d\n',X(n,1),X(n,2),s);
%
 %               break
  %          end
   %     end
   % end
%fclose(outputID);

[contour,N_s] = get_contour(n, DT, selected,adj_mat,invalid);

% %drop if the segmentation touches the boundary
% if length(unique([x1,x2])) > length(x1) 
%     x1 = [];
%     x2 = [];
%     y1 = [];
%     y2 = [];
%     return ; 
% end
% 
% if N_s ~= N_r && drop == true
%     disp('the simulation is dropped')
%     x1 = [];
%     x2 = [];
%     y1 = [];
%     y2 = [];
%     return ; 
% end


if draw_c == true    
        plot([x1';x2'], [y1';y2'],'Color',Color);
        axis([0 1 0 1])   
        axis square;
        hold on;
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

    function [sets_all_def, log_like_all_def] = merge_rand(num, sets_all_def, log_like_all_def,region_intensity,region_area,region_num_cells,cell_area, ...
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

        end
        
        
    end
end


