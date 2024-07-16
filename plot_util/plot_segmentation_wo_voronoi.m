
function [] = plot_segmentation_wo_voronoi(DT, selected, cx, cy, colors,pt_size,tri,seeds)

GRAY = [0.6 0.6 0.6];
hold on
if tri
    triplot(DT, 'Color', GRAY)
end
%hold on
% the final result
%selected = sets_all{index_BIC};
regionId = find(~cellfun(@isempty,selected));

index = 0;
for i = 1:length(regionId)
    index = index+1;
    scatter(cx(selected{regionId(i)}), cy(selected{regionId(i)}), pt_size,  colors(index, :), 'filled')
    %text(num2str(region_intensity(i)))
    axis equal
    axis([0 1 0 1])
end
%legend(string(regionId))
if nargin==8
    %plot seed
    for i = 1:length(seeds)
        scatter(cx(seeds{i}), cy(seeds{i}), 10,'K','filled')
    end
end
hold off
end



