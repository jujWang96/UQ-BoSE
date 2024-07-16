function [sets_all, log_like_all,region_intensity_all] = merge_region_random(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n,rand_num)
% region merging for an over-segmented graph randomly proportional to rank
% of BIC reduction. 
% Input variables:
%
% num: the number of regions in the over-segmented graph
% cell_area: the area of cells
% cell_log_intensity: the log intensity of cells
% region_sets: the region sets of the over-segmented graph
% adj_mat: the adjacent matrix of cells
% n: the number of cells
% rand_num: number of possible merge to choose each time
% m: penalty coefficient

% Output variables:
%
% sets_all: a series of region sets in the region merging process
% log_like_all: a series of log-likelihoods in the region merging process

% compute the intensity of regions and the connectivity among regions
[region_intensity, region_area, region_num_cells, adj_mat_region] =...
        get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);

    % compute the log-likelihood of the over-segmented graph before region merging
        log_like = get_metric_value('log-like', n, region_sets, region_area, region_intensity, region_num_cells);
        % construct a region merging tree
        % t represents the iter # of region merging
        region_indices = 1:num;
        log_like_all = zeros(num, 1);
        log_like_all(1) = log_like;
        sets_all = cell(num, 1);
        sets_all{1} = region_sets;
        for t = 1:num-1
            % find the pair of regions with the maximal log-likelihood increasement
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
            log_like_all(t+1) = log_like_all(t)+max_inc;
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
            tmp = sets_all{t};
            tmp{merge_region1} = [tmp{merge_region1} tmp{merge_region2}];
            tmp{merge_region2} = [];
            sets_all{t+1} = tmp;
            region_intensity_all{t+1} = region_intensity;

        end


end
