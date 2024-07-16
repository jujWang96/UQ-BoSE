function delta_log_like = get_delta_log_like(region_intensity_i, region_area_i,...
    region_num_cells_i, region_intensity_j, region_area_j, region_num_cells_j)
% Get the difference between the log-likelihood after merging region i and
% j and the log-likelihood before merging.

numerator = region_intensity_i*region_area_i+region_intensity_j*region_area_j;
sum_area = region_area_i+region_area_j;
delta_log_like = region_num_cells_i*log(numerator/region_intensity_i/sum_area)+...
    region_num_cells_j*log(numerator/region_intensity_j/sum_area);