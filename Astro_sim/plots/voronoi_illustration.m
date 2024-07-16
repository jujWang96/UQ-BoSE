rng(0)

n = 50;
X = rand([n, 2]);

h = figure;
cx = X(:, 1);
cy = X(:, 2);
DT = delaunayTriangulation(cx, cy);
triplot(DT, 'Color', 'k', 'LineWidth', 1.2)
hold on
scatter(cx, cy, 48, [0.8500, 0.3250, 0.0980], 'filled')
axis image
set(gca, 'fontsize', 16)
min_white_margin(gca, 0, 0, 0.02);
set(h, 'Position', [0, 0, 400, 350]);


h = figure;
voronoi(cx, cy);
set(findall(gca, 'Type', 'Line'), 'LineWidth', 1.2);
set(findall(gca, 'Type', 'Line'), 'Color', 'k');
hold on
scatter(cx, cy, 48, [0.8500, 0.3250, 0.0980], 'filled')
axis image
set(gca, 'fontsize', 16)
min_white_margin(gca, 0, 0, 0.02);
set(h, 'Position', [0, 0, 400, 350]);

