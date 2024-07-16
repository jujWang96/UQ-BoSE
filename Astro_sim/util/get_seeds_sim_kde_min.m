function [seeds, num_s] = get_seeds_sim_kde_min(cx, cy,invalid, adj_mat, size)
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
% size: number of seeds that will be initialized for each interval
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
step_size= length(idx_pdfX)/200;
% place seeds uniformly
seeds_all = [];
set_size=[];
for cdf_int = 1:step_size:length(idx_pdfX)-1

    if ceil(cdf_int+step_size)>max(idx_pdfX)
        photon_cdf_set = idx_pdfX(ceil(cdf_int):end);
    else
       photon_cdf_set = idx_pdfX(ceil(cdf_int):ceil(cdf_int+step_size)-1);
    end
       [~,index] = sort(pdfX(photon_cdf_set));
       if length(photon_cdf_set)<size
           seeds_all = [seeds_all, photon_cdf_set];
       else
           seeds_all = [seeds_all, photon_cdf_set(index(1:5))];
       end

end  



%reject if seed is in invalid 
seeds_all = setdiff(seeds_all,invalid);
seeds = num2cell(seeds_all);
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