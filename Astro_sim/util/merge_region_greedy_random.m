function [glb_selected,glb_region_intensity,glb_min_index,glb_min_BIC] = merge_region_greedy_random(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n,rand_step,rand_num,rep_num,m,reg_n)
%input: 
%  rand_num:number of candidate of each random merging step 
%  rand_step: number of steps of random merge
%  rep_num: number of repetition of random merge procedure 
%
%
%greedy merge the region until there are "rand_itr" number of region left 
%switch to random merge, each pick pick "rand_num" best candidate and
%randomly pick amoung them 
glb_min_BIC = realmax;
glb_min_index = 0;
%adjust the random merge step if there is not enough region 
if rand_step>=num-1
	rand_step = num-1;
end
[sets_all_def, log_like_all_def,region_intensity_def,region_area,region_num_cells,region_indices,adj_mat_region] = merge_def(num, cell_area, ...
        cell_log_intensity, region_sets, adj_mat, n,rand_step);
if rand_step==0
    BIC_all = -2*log_like_all_def+m*(num-1:-1:0)'*log(n);
    [glb_min_BIC, glb_min_index] = min(BIC_all);
     glb_selected = sets_all_def{glb_min_index};
    glb_region_intensity = region_intensity_def;
    return
end
for itr = 1:rep_num   
    [sets_all, log_like_all,region_intensity_all] = merge_rand(num, sets_all_def, log_like_all_def,region_intensity_def,region_area,region_num_cells, ...
     adj_mat_region,region_indices,rand_num,rand_step);
    BIC_all = -2*log_like_all+m*(num-1:-1:0)'*log(n);
    if nargin==10

    [min_BIC, min_index] = min(BIC_all);
    else
        min_index = num-reg_n+1;
        min_BIC = BIC_all(min_index);
    end
    if min_BIC<glb_min_BIC
        glb_min_BIC = min_BIC;
        glb_selected = sets_all{min_index};
        glb_region_intensity = region_intensity_all{min_index};
        glb_min_index = min_index;
    end
    
end

 function [sets_all_def, log_like_all_def,region_intensity,region_area,region_num_cells,region_indices,adj_mat_region] = merge_def(num, cell_area, ...
        cell_log_intensity, region_sets, adj_mat, n,rand_step)
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
        for t = 1:num-rand_step-1
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

    function [sets_all_def, log_like_all_def,region_intensity_all] = merge_rand(num, sets_all_def, log_like_all_def,region_intensity,region_area,region_num_cells, ...
     adj_mat_region,region_indices,rand_num,rand_itr)
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
