function [area_all, vx_all, vy_all] = get_area_vertex_for_tri_in_cell(edges, center, V)
% Get area and vertices for all the triangles in the Voronoi cell.
%
% Args:
%   edges: A M-by-2 matrix, each row of which is a pair of vertex indices
%     representing an edge of the Voronoi cell.
%   center: A 1-by-2 matrix representing the coordinate of the point that
%   expands the Voronoi cell. It will be paired with edge vertices to form
%     triangles.
%   V: A N-by-2 matrix, each row of which representing the coordinate of
%     a Voronoi cell vertex.
%
% Returns:
%   area_all: Area of all the triangles.
%   vx_all: num_tri-by-3 matrix, where each row gives the x-coordinates of
%     the three vertices of a triangle.
%   vy_all: num_tri-by-3 matrix, where each row gives the y-coordinates of
%     the three vertices of a triangle.

num_tri = size(edges, 1);
area_all = zeros(num_tri, 1);
vx_all = zeros(num_tri, 3);
vy_all = zeros(num_tri, 3);
for i = 1:num_tri
    v1 = V(edges(i, 1), :);
    v2 = V(edges(i, 2), :);
    v3 = center;
    vx_all(i, :) = [v1(1), v2(1), v3(1)];
    vy_all(i, :) = [v1(2), v2(2), v3(2)];
    area_all(i) = polyarea(vx_all(i, :), vy_all(i, :));
end

end