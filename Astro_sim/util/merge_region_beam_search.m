function [log_like_alls, sets_all,sets_opt, log_like_opt,region_area_bs,region_intensity_bs,region_num_cells_bs] = merge_region_beam_search(num, cell_area, ...
    cell_log_intensity, region_sets, adj_mat, n, bs,step)
% region merging for an over-segmented graph with beam search
% Input variables:
%
% num: the number of regions in the over-segmented graph
% cell_area: the area of cells
% cell_log_intensity: the log intensity of cells
% region_sets: the region sets of the over-segmented graph
% adj_mat: the adjacent matrix of cells
% n: the number of cells
% bs: beam search size 
% Output variables:
%
% sets_all: a series of region sets in the region merging process
% log_like_all: a series of log-likelihoods in the region merging process

% compute the intensity of regions and the connectivity among regions
[region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);

% compute the log-likelihood of the over-segmented graph before region merging
log_like = get_metric_value('log-like', n, region_sets, region_area, region_intensity, region_num_cells);
% construct region merging trees for beam search
% t represents the iter # of region merging
region_indices_bs = cell(bs,1);
for b = 1:bs
    region_indices_bs{b} = 1:num;
end
%seg_subset = repmat(mat2cell(1:num,[1],ones(1,num)),bs,1);
region_intensity_bs = repmat(region_intensity,1,bs);
region_num_cells_bs = repmat(region_num_cells,1,bs);
region_area_bs = repmat(region_area,1,bs);
adj_mat_region_bs = repmat(adj_mat_region,1,1,bs);
log_like_alls = zeros(num, bs);
log_like_alls(1,:) = log_like;
log_like_opt = zeros(num,1);
log_like_opt(1) = log_like;
sets_all = cell(num, bs);
sets_all(1,:) = {region_sets};
sets_opt = cell(num,1);
region_area_opt = cell(num,1);
region_intensity_opt = cell(num,1);
region_num_cells_opt = cell(num,1);
seg_subset = sets_all{1};
newseg_subset = sets_all{1};
for t = 1:step-1
    delta_log_like_candidate = -realmax*ones(bs,1);
    merge_candidate = cell(bs,1);
    tree_candidate = zeros(bs,1);
    currlog = -realmax*ones(bs,1);
    %start with a single segmentation
    if t == 1
        for i = region_indices_bs{1}
            for j = region_indices_bs{1}
                if i<j && adj_mat_region_bs(i, j,1)
                    numerator = region_intensity_bs(i,1)*region_area_bs(i,1)+region_intensity_bs(j,1)*region_area_bs(j,1);
                    sum_area = region_area_bs(i,1)+region_area_bs(j,1);
                    delta_log_like = region_num_cells_bs(i,1)*log(numerator/region_intensity_bs(i,1)/sum_area)+...
                        region_num_cells_bs(j,1)*log(numerator/region_intensity_bs(j,1)/sum_area);
                    % delta_log_like is defined as new log-like minus old
                    % log-like
                    % update the max delta value 
                    [max_inc, max_idx] = min(delta_log_like_candidate);
                    if delta_log_like>max_inc
                        %currseg = currseg_subset(1,:);
                        currseg  = seg_subset(1,:);
                        currseg{i} = [currseg{i},currseg{j}];
                        currseg{j} = [];
                        currlog(max_idx) = delta_log_like+log_like_alls(t,1);
                        newseg_subset(max_idx,:) = currseg;
                        delta_log_like_candidate(max_idx) = delta_log_like;
                        merge_candidate{max_idx} = [i,j];
                        tree_candidate(max_idx) = 1;
                    end
                end
            end
        end
    else
        %search for each candidate of bean search 
        for b = 1:bs
            for i = region_indices_bs{b}
                for j = region_indices_bs{b}
                    if i<j && adj_mat_region_bs(i, j,b)
                        numerator = region_intensity_bs(i,b)*region_area_bs(i,b)+region_intensity_bs(j,b)*region_area_bs(j,b);
                        sum_area = region_area_bs(i,b)+region_area_bs(j,b);
                        delta_log_like = region_num_cells_bs(i,b)*log(numerator/region_intensity_bs(i,b)/sum_area)+...
                            region_num_cells_bs(j,b)*log(numerator/region_intensity_bs(j,b)/sum_area);
                        % delta_log_like is defined as new log-like minus old
                        % log-like
                        % update the max delta value 
                        [min_max_log, min_idx] = min(currlog);
                        %
                        if delta_log_like+log_like_alls(t,b)>min_max_log

                            %check duplicity: 
                            duplicate = false;
                            currseg = seg_subset(b,:);
                            currseg{i} = [currseg{i},currseg{j}];
                            currseg{j} = [];
                            %check duplicity of loglikelihood and full
                            %segmentation
                            for c = 1:bs
                                if abs(currlog(c)-(delta_log_like+log_like_alls(t,b)) )<1e-10
                                    duplicate_seg = check_same_seg(newseg_subset(c,:),currseg,num);
                                   
                                    if duplicate_seg
                                        duplicate = true;
                                        break
                                    else
                                        continue
                                    end
                                end
                            end
                            if duplicate== false
                                currlog(min_idx) = delta_log_like+log_like_alls(t,b);
                                newseg_subset(min_idx,:) = currseg;
                                delta_log_like_candidate(min_idx) = delta_log_like;
                                merge_candidate{min_idx} = [i,j];
                                tree_candidate(min_idx) = b;
                            end
                        end
                    end
                end
            end
        end
    end
    %update trees

    %autofill the trees with the first tree 
    for b = 1:bs
        if tree_candidate(b)==0
            tree_candidate(b) = tree_candidate(1);
            newseg_subset(b,:) = newseg_subset(1,:);
            currlog(b)= currlog(1);
            merge_candidate{b} = merge_candidate{1};
            delta_log_like_candidate(b) = delta_log_like_candidate(1);
        end
    end
    
    region_indices_bs = region_indices_bs(tree_candidate);
    region_intensity_bs = region_intensity_bs(:,tree_candidate);
    region_num_cells_bs = region_num_cells_bs(:,tree_candidate);
    region_area_bs = region_area_bs(:,tree_candidate);
    adj_mat_region_bs = adj_mat_region_bs(:,:,tree_candidate);
    seg_subset = newseg_subset;
    
    max_log = -realmax;
    for b = 1:bs
        % update each subtree
        region_indices_bs{b}(region_indices_bs{b}==merge_candidate{b}(2)) = [];

        % update log-likelihood
        log_like_alls(t+1,b) = log_like_alls(t,tree_candidate(b))+delta_log_like_candidate(b);
  
        % merge merge_region1 and merge_region2
        % only keep merge_region1
        region_intensity_bs(merge_candidate{b}(1),b) = (region_intensity_bs(merge_candidate{b}(1),b)*region_area_bs(merge_candidate{b}(1),b)+...
            region_intensity_bs(merge_candidate{b}(2),b)*region_area_bs(merge_candidate{b}(2),b))/...
            (region_area_bs(merge_candidate{b}(1),b)+region_area_bs(merge_candidate{b}(2),b));
        region_area_bs(merge_candidate{b}(1),b) = region_area_bs(merge_candidate{b}(1),b)+region_area_bs(merge_candidate{b}(2),b);
        region_num_cells_bs(merge_candidate{b}(1),b) = region_num_cells_bs(merge_candidate{b}(1),b)+region_num_cells_bs(merge_candidate{b}(2),b);
        % update the connectivity
        index = setdiff(find(adj_mat_region_bs(merge_candidate{b}(2),:, b)==1), merge_candidate{b}(1));
        adj_mat_region_bs(merge_candidate{b}(1), index,b) = 1;
        adj_mat_region_bs(index, merge_candidate{b}(1),b) = 1;
        % record the sets after region merging
        tmp = sets_all{t,tree_candidate(b)};
        tmp{merge_candidate{b}(1)} = [tmp{merge_candidate{b}(1)} tmp{merge_candidate{b}(2)}];
        tmp{merge_candidate{b}(2)} = [];
        sets_all{t+1,b} = tmp;
        %store the best segmetation
        if log_like_alls(t+1,b)>max_log
            max_log = log_like_alls(t+1,b);
            log_like_opt(t+1) = max_log;
            %disp(log_like_opt(t+1))
            sets_opt{t+1} = sets_all{t+1,b};
            region_area_opt{t+1} = region_area_bs(:,b);
            region_intensity_opt{t+1} = region_intensity_bs(:,b);
            region_num_cells_opt{t+1} = region_num_cells_bs(:,b);
            sub_seg_opt = newseg_subset(b,:);
        end
    end
    

 
end

end