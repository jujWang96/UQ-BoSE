function [region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat)
% compute the intensity of regions and the connectivity among regions
%
% Input variables:
% num: number of regions
% cell_area: area of voronoi cells
% cell_log_intensity: log intensity of voronoi cells
% region_sets: (cell) region sets
% adj_mat: adjacent matrix of voronoi cells
%
% Output variables:
% region_intensity: intensity of regions
% region_area: area of regions
% region_num_cells: number of voronoi cells in each region
% adj_mat_region: adjacent matrix of regions


% compute the intensity
region_intensity = zeros(num, 1);
region_area = zeros(num, 1);
region_num_cells = zeros(num, 1);
for i = 1:num
    region_area(i) = sum(cell_area(region_sets{i}));
    region_intensity(i) = sum(cell_area(region_sets{i}).*exp(cell_log_intensity(region_sets{i})))/region_area(i);
    region_num_cells(i) = length(region_sets{i});
end

% compute the connectivity
adj_mat_region = zeros(num);
for i = 1:num-1
    for j = i+1:num
        adj_mat_region(i, j) = check_connect(i, j, region_sets, adj_mat);
        adj_mat_region(j, i) = adj_mat_region(i, j);
    end
end

end