function [sx_all, sy_all, acc_lambda] = get_bootstrap_samples(cx, cy,...
    DT, num_nonempty, selected_nonempty, log_int_all)
% Get bootstrap samples from fitted model.
% 
% Args:
%   cx: x-coordinates of photons.
%   cy: y-coordinates of photons.
%   DT: Delaunay triangulation object.
%   num_nonempty: Number of nonempty growing regions.
%   selected_nonempty: Matlab Cell of length num_nonempty. Each element is 
%     a list of Voronoi cell ids contained in the growing region.
%   log_int_all: Logarithm of intensity of all the growing regions.
%
% Returns:
%   sx_all: x-coordinates of bootstrap sampled points.
%   sy_all: y-coordinates of bootstrap sampled points.
%   acc_lambda: Accumulated (sum of) lambda.

[V, R] = voronoiDiagram(DT);
sx_all = [];
sy_all = [];
acc_lambda = 0;
for j = 1:num_nonempty
    cells_in_region = R(selected_nonempty{j});
    % Photons and Voronoi cells have the same ordering.
    cx_in_region = cx(selected_nonempty{j});
    cy_in_region = cy(selected_nonempty{j});
    intensity = exp(log_int_all(j));
    for i = 1:length(cells_in_region)
        % each row records two vertex indices of an edge
        edges = [cells_in_region{i}' [cells_in_region{i}(2:end)'; cells_in_region{i}(1)]];
        % sort vertex indices to avoid ambiguity
        edges = sort(edges, 2);
        center = [cx_in_region(i), cy_in_region(i)];
        [area_for_tri_in_cell, vx_for_tri_in_cell, vy_for_tri_in_cell] = get_area_vertex_for_tri_in_cell(edges, center, V);
        lambda_for_tri_in_cell = area_for_tri_in_cell * intensity;
        acc_lambda = acc_lambda + sum(lambda_for_tri_in_cell);
        num_for_tri_in_cell = poissrnd(lambda_for_tri_in_cell);
        % go through all the triangles in the cell
        for k = 1:length(num_for_tri_in_cell)
            % if there are points in the triangle
            if num_for_tri_in_cell(k) > 0
                vx = vx_for_tri_in_cell(k, :);
                vy = vy_for_tri_in_cell(k, :);
                [sx, sy] = uniform_sample_triangle(vx, vy, num_for_tri_in_cell(k));
                sx_all = [sx_all, sx];
                sy_all = [sy_all, sy];
            end
        end
    end
end

end