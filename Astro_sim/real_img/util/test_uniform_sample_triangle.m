vx = [0, 1, 0];
vy = [0, 0, 1];

[sx, sy] = uniform_sample_triangle(vx, vy, 1000);

scatter(sx, sy)
hold on
line([1, 0], [0, 1], 'Color', 'r')
axis equal
axis tight