function [] = plot_seeds2(DT, cx, cy, seeds, seeds_pt, seeds_rej, colors, num_s, num_s_pt)
%red 
GRAY = [0.6 0.6 0.6];

triplot(DT, 'Color', GRAY)
hold on
for i = 1:num_s
    scatter(cx(seeds{i}), cy(seeds{i}), 20, 'k', 's', 'filled')
end
if ~isempty(seeds_pt)
    for i = 1:num_s_pt
        scatter(cx(seeds_pt{i}), cy(seeds_pt{i}), 12, colors(1, :), 'd', 'filled')
    end
end
for i = 1:length(seeds_rej)
    scatter(cx(seeds_rej{i}), cy(seeds_rej{i}), 20, colors(2, :), '^', 'filled')
end
axis image

end