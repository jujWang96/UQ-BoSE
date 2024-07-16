function  [seeds, num_s] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat, cell_area,P,threshold)
% get initial seeds, which are either uniformly spreaded in the region or
% found by local maxima
%
% Input variables:
%
% 
% cx: the x coordinate of points
% cy: the y coordinate of points
% invalid: the set of points that are invalid
% adj_mat: the adjacent matrix. 1 means connected, 0 means not.
% cell_area: area of each voronoi cell
% P: number of percentiles 
% threshold: seeds with adjacent photon number <= threshold will not be initialized 
% Output variables:
%
% seeds: seed sets selected uniformly (excluding those rejected)
% seeds_rej: seed sets which are rejected
% seeds_pt: seed sets found by local maximum
% num_s: the number of seed sets selected uniformly (excluding those rejected)
% num_s_pt: the number of seed sets found by local maximum

n = length(cx);

% counter of number of seeds
num_s = 0;
% empty cell
seeds = [];

% set of points that have not been selected
unselected_points = 1:n;
[~, cell_sort_area] = sort(cell_area);


metric = 'l1';
curr_label = 1;
for k = 0:P-1
    %compute length scale LSCAL = 2*sqrt((area_q+area_(q+dq))/2)
    len_scale = 2*sqrt((prctile(cell_area,k/P*100)+prctile(cell_area,((k+1)/P)*100))/2);

    %disp(len_scale)
    %select all events such that {la_j are in percentile range [q,q+dq]}
    target_cell_set = cell_sort_area(max(ceil(k/P*n),1):ceil((k+1)/P*n));
    %initial the label = 0 for all cells 
    cell_label = zeros(1,length(target_cell_set));
    for idx = 1:length(target_cell_set)
        if cell_label(idx)>0
            continue
        else
            [min_label,candidate_in_range] = find_label_in_range(target_cell_set(idx),len_scale,target_cell_set,cell_label);
            %pause(0.5)
            if min_label==0
                cell_label(candidate_in_range) = curr_label;
                %scatter(cx(target_cell_set(candidate_in_range)),cy(target_cell_set(candidate_in_range)),12,colors(curr_label,:))

                curr_label = curr_label+1;
            else
                cell_label(candidate_in_range) = min_label;
                %scatter(cx(target_cell_set(candidate_in_range)),cy(target_cell_set(candidate_in_range)),12,colors(min_label,:))

            end

        end
        
    end
    for L = 1:curr_label
        candidate_cell_id = target_cell_set(cell_label==L);
        if length(candidate_cell_id)<threshold
            continue
        end
        num_s = num_s+1;
        seeds(num_s) = find_candidate_seed(candidate_cell_id,cell_area(candidate_cell_id),metric);
    end
    
    
    
end

%reject if seed is in invalid 
seeds = num2cell(setdiff(seeds,invalid));
num_s = length(seeds);
%seeds = num2cell(seeds_stable);
seeds_rej={};
seeds_pt = {};
num_s_pt = 0;


    function [min_label,candidate_in_range] = find_label_in_range(obs_idx,len_scale,target_cell_set,cell_label)
        %select all events within +-LSCAL of (x_j,y_j)
        
        minX = cx(obs_idx)-len_scale;
        maxX = cx(obs_idx)+len_scale;
        
        minY = cy(obs_idx)-len_scale;
        maxY = cy(obs_idx)+len_scale;
        %find smallest existing label LMIN
        candidate_in_range = cx(target_cell_set)>=minX & cx(target_cell_set)<=maxX ...
            & cy(target_cell_set)>=minY & cy(target_cell_set)<=maxY;
        assigned_label = cell_label(candidate_in_range);
        max_label = max(assigned_label);
        if max_label==0
            min_label = 0;
        else
            min_label = max_label;
            for l = assigned_label
                if l ~= 0 && l<min_label
                    min_label = l;
                end
            end
        end
    end

    function [min_cell_id]= find_candidate_seed(candidate_cell_id,candidate_cell_area,metric)
        label_set_size = length(candidate_cell_id);
        min_cell_id = candidate_cell_id(1);

        if label_set_size==1 
            return 
        elseif label_set_size==2
            [~,min_idx] = min(candidate_cell_area);
            min_cell_id = candidate_cell_id(min_idx);
            return
        else
            minval = Inf;
            %compute d_l = sum(distance from (x_l,y_l)) + area_l,
            %    select event kL which has smallest {d_l}

            D = pdist([cx(candidate_cell_id),cy(candidate_cell_id)],'euclidean');
            dist_mat = squareform(D);
            for cell_idx = 1:label_set_size
                currval = sum(dist_mat(cell_idx,:))+candidate_cell_area(cell_idx);
                if currval<minval
                    minval = currval;
                    min_cell_id = candidate_cell_id(cell_idx);
                end
                    
            end
                   
        end
    end

end
