function true_class = evaluate_points_sim4(points, loc, radius)

num_points = size(points, 1);
dist_mat = pdist2(points, loc);
true_class = zeros(num_points, 1);
for i = 1:num_points
    which_circle = find(dist_mat(i, :) <= radius);
    % The point is not in any of the circles.
    if isempty(which_circle)
        true_class(i) = 1;
    % The point is in the largest circle only.
    elseif length(which_circle) == 1
        true_class(i) = 2;
    % The point is in one of the smallest circles.
    else
        true_class(i) = sum(which_circle);
    end
end

end

