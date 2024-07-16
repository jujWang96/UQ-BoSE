function [cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, bound_x, bound_y, count)
% conduct some initial computations
% Input variables:
% 
% X: an n-by-2 matrix
% bound_x: the boundary of x-axis
% bound_y: the boundary of y-axis
% count: the count of photons at each location
%        for simulated data, it is one
%        for real data, it can be larger than one
%
% Output variables:
%
% cx: the first column of X
% cy: the second column of X
% n: the number of points
% DT: a Delaunay triangulation object
% E: a matrix of edges
% cell_log_intensity: the log intensity of voronoi cells
% cell_area: the area of voronoi cells

if nargin == 3
    n = size(X, 1);
    count = ones(n, 1);
end

cx = X(:, 1);
cy = X(:, 2);
DT = delaunayTriangulation(cx, cy);

[V, R] = voronoiDiagram(DT);
E = edges(DT);

n = length(cx);
intensity = zeros(n, 1);
cell_area = zeros(n, 1);
for i = 1:n
    % for the i-th polygon, obtain the max and min of x values and y values
    % for all vertices
    max_x = max(V(R{i}, 1));
    min_x = min(V(R{i}, 1));
    max_y = max(V(R{i}, 2));
    min_y = min(V(R{i}, 2));
    % check whether any of the vertices is out of the boundary
    % if so, assign NaN to cell_area since the estimate is not reliable
    if max_x>bound_x(2) || min_x<bound_x(1) || max_y>bound_y(2) || min_y<bound_y(1)
        cell_area(i) = NaN;
        intensity(i) = NaN;
    else
        cell_area(i) = polyarea(V(R{i}, 1), V(R{i}, 2));
        intensity(i) = count(i)/cell_area(i);
    end
end

cell_log_intensity = log(intensity);

end