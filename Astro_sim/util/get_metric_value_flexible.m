function value = get_metric_value_flexible(metric, N, region_sets, region_area, region_intensity, region_num_cells)
% compute the value of metric allowing empty array in region_set
%
% Input variables:
%
% metric: metric of interest
% N: total number of points
% region_sets: region sets
% region_area: area of regions
% region_intensity: intensity of regions
% region_num_cells: number of cells in each region
%
% Output variables:
%
% value: value of the metric

% number of regions
K = length(region_sets);
%log(n!)
log_like = -sum(log(1:N));
nonempty = find(~cellfun(@isempty,region_sets));
for k = nonempty
    
    log_like = log_like+region_num_cells(k)*log(region_intensity(k))-...
        region_intensity(k)*region_area(k);
end


if strcmp(metric, 'BIC')
    % use number of regions as number of free parameters
    % this model complexity penalty is not large enough
    value = -2*log_like+K*log(N);
elseif strcmp(metric, 'log-like')
    value = log_like;
end
    
end