function [] = plot_circles(loc, radius,lw,colr)

if nargin==3
    colr = [0.6 0.6 0.6];
end
for i = 1:length(radius)
    viscircles(loc(i, :), radius(i), 'EdgeColor', colr, 'LineWidth', lw);
end

end
