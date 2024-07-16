function [] = plot_segmentation_select(DT, selected, cx, cy, colors)

GRAY = [0.6 0.6 0.6];
%figure
triplot(DT, 'Color', GRAY)
hold on
% the final result
index = 0;
for i = 1:length(selected)
    if ~isempty(selected{i})
        index = index+1;
        scatter(cx(selected{i}), cy(selected{i}),8, colors(index, :),'.');
    end
        axis equal

    axis([0 1 0 1])
end

end
