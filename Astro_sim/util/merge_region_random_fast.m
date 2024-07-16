function [sets_all, log_like_all] = merge_region_random_fast(num_region, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n, rand_num, seed)
% Fast region merging with randomization for an over-segmented graph.
%
% Args:
%    num_region: Number of regions in the over-segmented graph.
%    region_area: Area of all regions.
%    region_intensity: Intensity of all regions.
%    region_sets: Region sets of the over-segmented graph.
%    adj_mat_region: Adjacent matrix of regions.
%    region_num_cells: Number of cells in each region.
%    n: Number of data points.
%    rand_num: Number of candidate to randomly selected from for each state
%
% Returns:
%   sets_all: A series of region sets during the region merging process.
%   log_like_all: A series of log-likelihoods during the region merging 
%       process.

% Compute the log-likelihood of the over-segmented graph before region
% merging.
if nargin == 9
    rng(seed)
end
log_like = get_metric_value('log-like', n, region_sets, region_area, region_intensity, region_num_cells);

log_like_all = zeros(num_region, 1);
log_like_all(1) = log_like;
sets_all = cell(num_region, 1);
sets_all{1} = region_sets;

% Get all adjacent region pairs.
[region_ind1, region_ind2] = find(adj_mat_region);
% Each row of delta_log_like records the change of the log-likelihood
% (first element) if merging the two neighboring regions (second and third
% elements).
% Due to the symmetry of the adjacent matrix, only use half of all the
% pairs.
delta_log_like = zeros(length(region_ind1)/2, 3);
index = 0;
for k = 1:length(region_ind1)
    i = region_ind1(k);
    j = region_ind2(k);
    if i < j
        index = index + 1;
        delta_log_like(index, 3) = get_delta_log_like(region_intensity(i), region_area(i),...
            region_num_cells(i), region_intensity(j), region_area(j), region_num_cells(j));
        delta_log_like(index, 1) = i;
        delta_log_like(index, 2) = j;
    end
end

for t = 1:num_region-1
%     if mod(t, 100) == 0 
%         disp(t)
%         disp(size(delta_log_like, 1))
%     end
    % Find the pair of regions with the top "num_rand" log-likelihood
    % increasement.
    max_inc_candidate = -realmax*ones(1,rand_num);
    merge_region1_candidate = zeros(1,rand_num);
    merge_region2_candidate = zeros(1,rand_num);
        
    for k = 1:size(delta_log_like, 1)
        i = delta_log_like(k, 1);
        j = delta_log_like(k, 2);
        current_delta_log_like = delta_log_like(k, 3);
        % The larger, the better.
        [min_max_inc, min_idx] = min(max_inc_candidate);
        if current_delta_log_like>min_max_inc
            max_inc_candidate(min_idx) = current_delta_log_like;
            merge_region1_candidate(min_idx) = i;
            merge_region2_candidate(min_idx) = j;
        end
    end
    %randomly select one valid candidate to merge (edge case: possible segmentation < rand_num)
   
    if size(delta_log_like, 1)>0 
        %ensure that there is at least one candidate left to choose from 
        valid_merge = max_inc_candidate>-realmax;
        max_inc_candidate = max_inc_candidate(valid_merge);
        merge_region1_candidate = merge_region1_candidate(valid_merge);
        merge_region2_candidate = merge_region2_candidate(valid_merge);
        idx = randi([1,min(rand_num,length(max_inc_candidate))],1);
        max_inc = max_inc_candidate(idx);
    
        merge_region1 = merge_region1_candidate(idx);
        merge_region2 = merge_region2_candidate(idx);
    end
    % Any rows containing merge_region1 or merge_region2 are not valid.
    rows_to_delete = delta_log_like(:, 1) == merge_region1 |...
        delta_log_like(:, 2) == merge_region1 |...
        delta_log_like(:, 1) == merge_region2 |...
        delta_log_like(:, 2) == merge_region2;
    delta_log_like(rows_to_delete, :) = [];
    % Update log-likelihood.
    log_like_all(t+1) = log_like_all(t)+max_inc;
    % Merge merge_region1 and merge_region2.
    % Only keep merge_region1.
    region_intensity(merge_region1) = (region_intensity(merge_region1)*region_area(merge_region1)+...
        region_intensity(merge_region2)*region_area(merge_region2))/...
        (region_area(merge_region1)+region_area(merge_region2));
    region_area(merge_region1) = region_area(merge_region1)+region_area(merge_region2);
    region_num_cells(merge_region1) = region_num_cells(merge_region1)+region_num_cells(merge_region2);
    % Update the connectivity.
    % Find all neighbors of merge_region2.
    index = setdiff(find(adj_mat_region(merge_region2, :)==1), merge_region1);
    % They become neighbors of merge_region1 after merging.
    adj_mat_region(merge_region1, index) = 1;
    adj_mat_region(index, merge_region1) = 1;
    % Deactivate merge_region2.
    adj_mat_region(merge_region2, :) = 0;
    adj_mat_region(:, merge_region2) = 0;
    % Update delta_log_like. Note that any rows containing merge_region1 or 
    % merge_region2 have been deleted. So we need to add back rows 
    % containing merge_region1.
    % Find all neighbors of merge_region1 after updating EVERYTHING.
    region_ind2 = find(adj_mat_region(merge_region1, :));
    add_delta_log_like = zeros(length(region_ind2), 3);
    for j = 1:length(region_ind2)
        add_delta_log_like(j, 1) = merge_region1;
        add_delta_log_like(j, 2) = region_ind2(j);
        add_delta_log_like(j, 3) = get_delta_log_like(region_intensity(merge_region1), region_area(merge_region1),...
            region_num_cells(merge_region1), region_intensity(region_ind2(j)), region_area(region_ind2(j)),...
            region_num_cells(region_ind2(j)));
    end
    delta_log_like = [delta_log_like; add_delta_log_like];
    % Record the sets after region merging.
    tmp = sets_all{t};
    tmp{merge_region1} = [tmp{merge_region1} tmp{merge_region2}];
    tmp{merge_region2} = [];
    sets_all{t+1} = tmp;
end

end