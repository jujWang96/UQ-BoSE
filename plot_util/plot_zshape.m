function [] = plot_zshape(range_x, range_y,lw,colr)

if nargin == 3
    colr = [0.6 0.6 0.6];
end

% horizontal lines
for i = 1:3
    plot(range_x(i, :), [range_y(i, 1) range_y(i, 1)], 'Color', colr, 'LineWidth', lw)
    hold on;
    plot(range_x(i, :), [range_y(i, 2) range_y(i, 2)], 'Color', colr, 'LineWidth', lw)
    hold on;
end
% vertical lines
plot([range_x(1, 1) range_x(1, 1)], range_y(1, :), 'Color', colr, 'LineWidth', lw)
hold on;
plot([range_x(3, 2) range_x(3, 2)], range_y(3, :), 'Color', colr, 'LineWidth', lw)
hold on;

plot([range_x(1, 2) range_x(1, 2)], [range_y(1, 2) range_y(3, 2)], 'Color', colr, 'LineWidth', lw)
hold on;

plot([range_x(2, 2) range_x(2, 2)], [range_y(1, 1) range_y(3, 1)], 'Color', colr, 'LineWidth', lw)
hold on;


end