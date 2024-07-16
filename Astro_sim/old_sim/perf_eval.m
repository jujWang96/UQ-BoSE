function [dr, far, err, region_area, source_x, source_y, n_S_correct] =...
    perf_eval(n_source, loc, radius, cx, cy, selected, cell_log_intensity, cell_area)
% performance evaluation

% compute/estimate the intensity of the fitted regions
region_area = zeros(n_source+1, 1);
region_intensity = zeros(n_source+1, 1);
for i = 1:n_source+1
    region_area(i) = sum(cell_area(selected{i}));
    region_intensity(i) = sum(cell_area(selected{i}).*exp(cell_log_intensity(selected{i})))/region_area(i);
end
[background_intensity, index_background] = min(region_intensity);
% remove background
selected(index_background) = [];
region_area(index_background) = [];
% n1 records the number of successfully detected source photons for
% each source
% n2 records the number of false alarms for each source
n1 = zeros(n_source, 1);
n2 = zeros(n_source, 1);
% n_S records the number of photons in S_i
% n_Sc records the number of photons in S_i^c
% also include photons from the background
n_S = zeros(n_source, 1);
n_Sc = zeros(n_source, 1); 
% estimated source position
source_x = zeros(n_source, 1);
source_y = zeros(n_source, 1);
n_S_correct = zeros(n_source, 1);
% note that there might be a mismatch between estimated and true sources
for i = 1:n_source 
    dist_sq = (cx-loc(i, 1)).^2+(cy-loc(i, 2)).^2;
    n_S(i) = length(find(dist_sq<=radius(i)^2));
    n_Sc(i) = length(find(dist_sq>radius(i)^2));
    denom = sum(exp(cell_log_intensity(selected{i})));
    source_x(i) = sum(exp(cell_log_intensity(selected{i})).*cx(selected{i}))/denom;
    source_y(i) = sum(exp(cell_log_intensity(selected{i})).*cy(selected{i}))/denom;
    % with background correction
    % min(region_intensity) is supposed to be the intensity of the
    % background
    n_S_correct(i) = length(selected{i})-background_intensity*region_area(i);
end
index_set = 1:n_source;
% record how estimated and true sources match with each other
match = zeros(n_source, 1);
for i = 1:n_source
    % for each estimated source, pick the closest true source to match
    dist_sq = (source_x(i)-loc(index_set, 1)).^2+(source_y(i)-loc(index_set, 2)).^2;
    [~, index_tmp] = min(dist_sq);
    % the index of matched true source
    index_true_source = index_set(index_tmp);
    % delete it from index_set
    index_set(index_tmp) = [];
    % record it
    match(i) = index_true_source;
end

% compute n1 and n2 after matching
for i = 1:n_source
    dist_sq = (cx(selected{i})-loc(match(i), 1)).^2+(cy(selected{i})-loc(match(i), 2)).^2;
    n1(i) = length(find(dist_sq<=radius(match(i))^2));
    n2(i) = length(find(dist_sq>radius(match(i))^2));
end
  
% detection rate
dr = n1./n_S(match);
% false alarm rate
far = n2./n_Sc(match);
% total misclassification rate
err = 1-dr+far;
   
end