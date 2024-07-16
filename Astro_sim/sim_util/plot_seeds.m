function [] = plot_seeds(DT, cx, cy, seeds, seeds_pt, seeds_rej, colors, num_s, num_s_pt)

GRAY = [0.6 0.6 0.6];

triplot(DT, 'Color', GRAY)
hold on
for i = 1:num_s
    scatter(cx(seeds{i}), cy(seeds{i}), 12, colors(i, :), 's', 'filled')
end
if ~isempty(seeds_pt)
    for i = 1:num_s_pt
        scatter(cx(seeds_pt{i}), cy(seeds_pt{i}), 12, colors(i + num_s, :), 'd', 'filled')
    end
end
for i = 1:length(seeds_rej)
    scatter(cx(seeds_rej{i}), cy(seeds_rej{i}), 12, 'k', '^', 'filled')
end
axis image

end