% test function for evaluate_points_sim6
num_points = 10000;
rng(0)
cx = rand(num_points, 1);
cy = rand(num_points, 1);

loc_ring = [0.3 0.5];
radius_out = 0.4;
radius_in = 0.2;

loc = [0.5 0.7; 0.6 0.5; 0.5 0.3];
radius = 0.025*ones(1, 3);
true_class = evaluate_points_sim6([cx cy], loc_ring, radius_out, radius_in, loc, radius);

scatter(cx, cy, 12, true_class, 'filled')
colorbar
colormap(jet)
axis image