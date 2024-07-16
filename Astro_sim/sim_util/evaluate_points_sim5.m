function true_class = evaluate_points_sim5(points, range_x, range_y, loc, radius)

num_points = size(points, 1);
num_rectangles = size(range_x, 1);
% i-th row and j-th col means whether the i-th point is in the j-th
% rectangle
in_all = zeros(num_points, num_rectangles);
for i = 1:num_rectangles
    in_all(:, i) = all([points(:, 1) >= range_x(i, 1) points(:, 1) <= range_x(i, 2)...
        points(:, 2) >= range_y(i, 1) points(:, 2) <= range_y(i, 2)], 2);
end
% whether the i-th point is in any of the rectangles
in = any(in_all, 2);
dist_mat = pdist2(points, loc);
true_class = zeros(num_points, 1);
for i = 1:num_points
    % If it is not empty, it would have only one element.
    which_circle = find(dist_mat(i, :) <= radius, 1);
    % The point is not in any of the circles.
    if isempty(which_circle)
        % in the L shape
        if in(i)
            true_class(i) = 2;
        else
            true_class(i) = 1;
        end
    else
        % starts from 3
        true_class(i) = which_circle + 2;
    end
end

end

