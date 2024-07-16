function [] = plot_seeds_one(DT, cx, cy, seeds, seeds_pt, seeds_rej, colors, num_s, num_s_pt)
%plot seeds by one seed in the middle of the cluster 
GRAY = [0.6 0.6 0.6];

triplot(DT, 'Color', GRAY)
hold on
for i = 1:num_s
    [cxc,cyc] = find_centroid(cx(seeds{i}),cy(seeds{i}));
    scatter(cxc, cyc, 30, 'k', 's', 'filled')

end
if ~isempty(seeds_pt)
    for i = 1:num_s_pt
        [cxc,cyc] = find_centroid(cx(seeds_pt{i}),cy(seeds_pt{i}));

        scatter(cxc,cyc, 30, colors(1, :), 'd', 'filled')
    end
end
for i = 1:length(seeds_rej)
    [cxc,cyc] = find_centroid(cx(seeds_rej{i}),cy(seeds_rej{i}));

    scatter(cxc, cyc, 30, colors(2, :), '^', 'filled')
end
axis image



    %find the centroid photon of cx and cy
    function [cxc,cyc] = find_centroid(cx,cy)
        dxv = abs(bsxfun(@minus,cx,cx'));
        dyv = abs(bsxfun(@minus,cy,cy'));
        dsum = sum(dxv+dyv);
        [~,I] = min(dsum);
        cxc = cx(I);
        cyc = cy(I);
    end


end