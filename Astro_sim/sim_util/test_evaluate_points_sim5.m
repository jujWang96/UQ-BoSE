% test function for evaluate_points_sim5
num_points = 10000;
rng(0)
cx = rand(num_points, 1);
cy = rand(num_points, 1);

range_x = [0.2 0.4; 0.4 0.6; 0.6 0.8];
range_y = [0.2 0.4; 0.2 0.8; 0.6 0.8];
loc = [0.35 0.3; 0.5 0.5; 0.65 0.7];
radius = [0.025 0.025 0.025];
true_class = evaluate_points_sim5([cx cy], range_x, range_y, loc, radius);

scatter(cx, cy, 12, true_class, 'filled')
colorbar
colormap(jet)
axis image