function true_class = evaluate_points_sim6(points, loc_ring, radius_out, radius_in, loc, radius)

num_points = size(points, 1);
dist_mat_ring = pdist2(points, loc_ring);
dist_mat = pdist2(points, loc);
true_class = zeros(num_points, 1);
for i = 1:num_points
    in_ring = (dist_mat_ring(i) >= radius_in) & (dist_mat_ring(i) <= radius_out) &...
        (points(i, 1) >=loc_ring(1));
    if ~in_ring
        true_class(i) = 1;
    else
        which_circle = find(dist_mat(i, :) <= radius, 1);
        if isempty(which_circle)
            true_class(i) = 2;
        else
            true_class(i) = which_circle + 2;
        end
    end        
end

end

