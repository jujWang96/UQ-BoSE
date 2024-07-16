function [] = plot_segmentation(DT, index_BIC, sets_all, cx, cy, colors)

GRAY = [0.6 0.6 0.6];

triplot(DT, 'Color', GRAY)
hold on
% the final result
selected = sets_all{index_BIC};
index = 0;
for i = 1:length(selected)
    if ~isempty(selected{i})
        index = index+1;
        scatter(cx(selected{i}), cy(selected{i}), 12,  colors(index, :), 'filled')
    end
end

end
