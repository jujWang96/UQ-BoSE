% test function for evaluate_points_sim4
num_points = 10000;
rng(0)
cx = rand(num_points, 1);
cy = rand(num_points, 1);

loc = [0.5 0.5; 0.4 0.4; 0.4 0.6; 0.6 0.4; 0.6 0.6];
radius = [0.25 0.025*ones(1, 4)];
true_class = evaluate_points_sim4([cx cy], loc, radius);

scatter(cx, cy, 12, true_class, 'filled')
colorbar
colormap(jet)
axis image