function [log_like_buffer, sets_buffer, sets_opt, log_like_opt, region_area_buffer,region_intensity_buffer, region_num_cells_buffer] = merge_region_fast_beamsearch_optimization(num_region, region_area, ...
    region_intensity, region_sets, adj_mat_region, region_num_cells, n, bw)
% Fast region merging for an over-segmented graph with beam merge.
%
% Args:
%    num_region: Number of regions in the over-segmented graph.
%    region_area: Area of all regions.
%    region_intensity: Intensity of all regions.
%    region_sets: Region sets of the over-segmented graph.
%    adj_mat_region: Adjacent matrix of regions.
%    region_num_cells: Number of cells in each region.
%    n: Number of data points.
%    bs: beam width of the buffer 
% Returns:
%   sets_opt: A series of optimum region sets during the region merging process.
%   log_like_opt: A series of optimum log-likelihoods during the region merging 
%       process.
%   log_like_buffer:A series of buffer of log-likelihoods for each state
%   sets_buffer: A series of buffer of region sets for each state
%   

% Compute the log-likelihood of the over-segmented graph before region
% merging.
log_like = get_metric_value('log-like', n, region_sets, region_area, region_intensity, region_num_cells);

log_like_opt = zeros(num_region, 1);
log_like_opt(1) = log_like;
sets_opt = cell(num_region, 1);
sets_opt{1} = region_sets;
region_intensity_buffer = repmat(region_intensity,1,bw);
region_num_cells_buffer = repmat(region_num_cells,1,bw);
region_area_buffer = repmat(region_area,1,bw);
adj_mat_region_buffer = repmat(adj_mat_region,1,1,bw);
sets_buffer = cell(num_region, bw);
sets_buffer(1,:) = {region_sets};
log_like_buffer = -realmax*ones(num_region, bw);
log_like_buffer(1,1) = log_like;

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

delta_log_like_buffer = cell(1,1);
delta_log_like_buffer(:) = {delta_log_like};

curr_bs = 1;
for t = 1:num_region-1
%     if mod(t, 100) == 0 
%         disp(t)
%         disp(size(delta_log_like, 1))
%     end
  
    %for each buffer search for the top bs pairs of candidate
    
    max_log_candidate = -realmax*ones(1,bw);
        merge_region1_candidate = zeros(1,bw);
        merge_region2_candidate = zeros(1,bw);
        parent_buffer_candidate = zeros(1,bw);
        seg_cell_buffer = cell(bw,num_region);
    for b = 1:curr_bs
        
        for k = 1:size(delta_log_like_buffer{b}, 1)
            i = delta_log_like_buffer{b}(k, 1);
            j = delta_log_like_buffer{b}(k, 2);
            
            current_delta_log_like = delta_log_like_buffer{b}(k, 3);
            % The larger, the better.
             [max_log, max_idx] = min(max_log_candidate);
             curr_log = current_delta_log_like+log_like_buffer(t,b);
            if curr_log>max_log
                %check duplicate segmentation 
                tmp = sets_buffer{t,b};
                tmp{i} = sort([tmp{i} tmp{j}]);
                tmp{j} = [];
                %tmp_segstr = seg_to_str(tmp,num_region);
                duplicate = false;
                for c = 1:bw
                    if abs(max_log_candidate(c)-curr_log)<1e-10
%                         if 
%                             duplicate = true;
%                             break;

                        %debug 
%                         if (strcmp(tmp_segstr,seg_curr_buffer(c)) && ~check_same_seg(seg_cell_buffer(c,:),tmp,num_region))...
%                             || (~strcmp(tmp_segstr,seg_curr_buffer(c)) && check_same_seg(seg_cell_buffer(c,:),tmp,num_region))
%                                 disp('err')
%                         end
                        %if strcmp(tmp_segstr,seg_curr_buffer(c))
                        if check_same_seg_restrict(seg_cell_buffer(c,:),tmp)    
                            duplicate = true;
                            break
                        else
                            continue
                        end
                    end
                end
                
                if ~duplicate
                    %seg_curr_buffer(max_idx) = tmp_segstr;
                    seg_cell_buffer(max_idx,:) = tmp;
                    max_log_candidate(max_idx) = curr_log;
                    merge_region1_candidate(max_idx) = i;
                    merge_region2_candidate(max_idx) = j;
                    parent_buffer_candidate(max_idx) = b;
                 %disp(buffer_set)
%                     if ~ismember(tmp_segstr,buffer_set)
%                         buffer_set(max_idx) = tmp_segstr;
%                         
%                     end 
                end
            end
        end
      
    end
    
    %select the unique segmentation with top bs loglikelihood 
    %drop if there is -realmax value contained (buffer size < bs)

    idx_buffer = find(max_log_candidate~=-realmax);
    %select the top bs non repeated segmentations
    merge_region1 = reshape(merge_region1_candidate(idx_buffer),1,[]);
    merge_region2 = reshape(merge_region2_candidate(idx_buffer),1,[]);
    parent_buffer = reshape(parent_buffer_candidate(idx_buffer),1,[]);

    max_log_buffer = reshape(max_log_candidate(idx_buffer),1,[]);

    curr_bs = min(length(idx_buffer),bw);

    % Merge merge_region1 and merge_region2 for each element.
    % Only keep merge_region1.
    
    region_intensity_buffer_new = region_intensity_buffer(:,parent_buffer);
    region_area_buffer_new = region_area_buffer(:,parent_buffer);
    region_num_cells_buffer_new = region_num_cells_buffer(:,parent_buffer);
       region_intensity_buffer_new(merge_region1+(0:curr_bs-1)*num_region) = (region_intensity_buffer_new(merge_region1+(0:curr_bs-1)*num_region).*region_area_buffer_new(merge_region1+(0:curr_bs-1)*num_region)+...
        region_intensity_buffer_new(merge_region2+(0:curr_bs-1)*num_region).*region_area_buffer_new(merge_region2+(0:curr_bs-1)*num_region))./...
        (region_area_buffer_new(merge_region1+(0:curr_bs-1)*num_region)+region_area_buffer_new(merge_region2+(0:curr_bs-1)*num_region));
    region_area_buffer_new(merge_region1+(0:curr_bs-1)*num_region) = region_area_buffer_new(merge_region1+(0:curr_bs-1)*num_region)+region_area_buffer_new(merge_region2+(0:curr_bs-1)*num_region);
    region_num_cells_buffer_new(merge_region1+(0:curr_bs-1)*num_region) = region_num_cells_buffer_new(merge_region1+(0:curr_bs-1)*num_region)+region_num_cells_buffer_new(merge_region2+(0:curr_bs-1)*num_region);


    
    % Update the connectivity.
    % Find all neighbors of merge_region2.
    delta_log_like_buffer_new = cell(1,curr_bs);
    adj_mat_region_buffer_new = repmat(adj_mat_region,1,1,curr_bs);
     for k = 1:curr_bs
         % Any rows containing merge_region1 or merge_region2 are not valid.
        b = parent_buffer(k);
        rows_to_delete = delta_log_like_buffer{b}(:, 1) == merge_region1(k) |...
        delta_log_like_buffer{b}(:, 2) == merge_region1(k) |...
        delta_log_like_buffer{b}(:, 1) == merge_region2(k) |...
        delta_log_like_buffer{b}(:, 2) == merge_region2(k);
        delta_log_like_buffer_new{k} = delta_log_like_buffer{b};
      
        delta_log_like_buffer_new{k}(rows_to_delete, :) = [];

        index = setdiff(find(adj_mat_region_buffer(merge_region2(k), :,b)==1), merge_region1(k));
        % They become neighbors of merge_region1 after merging.
        adj_mat_region_buffer_new(:,:,k) = adj_mat_region_buffer(:,:,b);
        adj_mat_region_buffer_new(merge_region1(k), index,k) = 1;
        adj_mat_region_buffer_new(index, merge_region1(k),k) = 1;
        % Deactivate merge_region2.
        adj_mat_region_buffer_new(merge_region2(k), :,k) = 0;
        adj_mat_region_buffer_new(:, merge_region2(k),k) = 0;
        % Update delta_log_like. Note that any rows containing merge_region1 or 
        % merge_region2 have been deleted. So we need to add back rows 
        % containing merge_region1.
        % Find all neighbors of merge_region1 after updating EVERYTHING.
        % Update adjacent matrix
        region_ind2 = find(adj_mat_region_buffer_new(merge_region1(k), :, k));
        add_delta_log_like = zeros(length(region_ind2), 3);
        for j = 1:length(region_ind2)
            %always merge larger to smaller index!
            add_delta_log_like(j, 1) = min(merge_region1(k),region_ind2(j));
            add_delta_log_like(j, 2) = max(merge_region1(k),region_ind2(j));
            add_delta_log_like(j, 3) = get_delta_log_like(region_intensity_buffer_new(merge_region1(k),k), region_area_buffer_new(merge_region1(k),k),...
                region_num_cells_buffer_new(merge_region1(k),k), region_intensity_buffer_new(region_ind2(j),k), region_area_buffer_new(region_ind2(j),k),...
                region_num_cells_buffer_new(region_ind2(j),k));
        end
        delta_log_like_buffer_new{k} = [delta_log_like_buffer_new{k}; add_delta_log_like];
        % Record the sets after region merging.
%         tmp = sets_buffer{t,b};
%         tmp{merge_region1(k)} = sort([tmp{merge_region1(k)} tmp{merge_region2(k)}]);
%         tmp{merge_region2(k)} = [];
%         sets_buffer{t+1,k} = tmp;

     end
     for b = 1:curr_bs
        sets_buffer{t+1,b} = seg_cell_buffer(b,:);
     end
    % Update log-likelihood for each element in the buffer 
    delta_log_like_buffer = delta_log_like_buffer_new;
    adj_mat_region_buffer = adj_mat_region_buffer_new;
    log_like_buffer(t+1,1:length(max_log_buffer)) = max_log_buffer;

    [max_log,max_idx] = max(log_like_buffer(t+1,:));
        
    log_like_opt(t+1) = max_log;
    sets_opt{t+1} = sets_buffer{t+1,max_idx};
    region_intensity_buffer = region_intensity_buffer_new;
    region_area_buffer = region_area_buffer_new;
    region_num_cells_buffer = region_num_cells_buffer_new;
    %check 

end
end
