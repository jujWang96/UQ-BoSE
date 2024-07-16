function [seeds, num_s] = get_seeds_sim_kde(cx, cy,invalid, adj_mat, P,threshold)
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
seeds = {};

% set of points that have not been selected
unselected_points = 1:n;
unselected_points = setdiff(unselected_points, invalid);
cx_valid = cx(unselected_points);
cy_valid = cy(unselected_points);
%place seeds based on the each interval of density 
%[cdfX, Xi]= ksdensity([cx,cy],[cx_valid,cy_valid]);
[pdfX, Xi] = ksdensity([cx,cy],[cx,cy]);
[~,idx_pdfX] = sort(pdfX);
step_size= length(idx_pdfX)/P;
% place seeds uniformly
seeds_all = [];
set_size=[];
for cdf_int = 1:step_size:length(idx_pdfX)-1

    if ceil(cdf_int+step_size)>max(idx_pdfX)
        photon_cdf_set = idx_pdfX(ceil(cdf_int):end);
    else
       photon_cdf_set = idx_pdfX(ceil(cdf_int):ceil(cdf_int+step_size)-1);
    end
       adj_mat_seed = adj_mat(photon_cdf_set, photon_cdf_set);
       %record the slices of the inte
        bins = conncomp(graph(adj_mat_seed));
        for bin = unique(bins)
            %initialize one seed on the photon nearest to the centroid of
            %the connected region.
            candidate_seed_set = photon_cdf_set(bins==bin);
            %scatter(cx(candidate_seed_set),cy(candidate_seed_set),'.')

            set_center = mean([cx(candidate_seed_set),cy(candidate_seed_set)],1);
            dist = (cx(candidate_seed_set)-set_center(1)).^2+...
                (cy(candidate_seed_set)-set_center(2)).^2;
            [~, index] = sort(dist);
            seeds_all = [seeds_all candidate_seed_set(index(1))];
            scatter(cx(candidate_seed_set(index(1))),cy(candidate_seed_set(index(1))))
            set_size = [set_size sum(bins==bin)];
        end 
        axis equal
        axis([0 1 0 1])
end  


%initialize seeds only if adjacent photon number larger than threshold 
seeds_stable = seeds_all(set_size>threshold);
%reject if seed is in invalid 
seeds_stable = setdiff(seeds_stable,invalid);
seeds = num2cell(seeds_stable);
num_s = length(seeds);

% % seed rejection based on area after initial selection
% index = [];
% for i = 1:num_s
%     areas = cell_area(seeds{i});
%     lambda_inv = mean(areas);
%     std_area = 0.53*lambda_inv;
%     max_range = 2*factor*std_area;
%     seed_range = max(areas)-min(areas);
%     if seed_range>max_range
%         index = [index i];
%     end
% end
% % remove
% seeds_rej = seeds(index);
% seeds(index) = [];
% num_s = num_s-length(index);


end