real_full_fast
figure
%triplot(DT, 'Color', GRAY)
hold on
init_region = sets_all{index_BIC};
lbd = 0.0003;
for i = 1:num
    if ~isempty(init_region{i})
      log_int = log(sum(exp(cell_log_intensity(init_region{i})).*cell_area(init_region{i}))/sum(cell_area(init_region{i})));
       if sum(cell_area(init_region{i}))<lbd
            scatter(cx(init_region{i}), cy(init_region{i}), 12, log_int*ones(length(init_region{i}), 1), 'filled')
        end
    end
end
cb = colorbar('EastOutside');
%colormap(hsv);
colormap default
%colormap(turbo)
%colormap(jet)

cb.Ticks = linspace(12, 18, 7);

axis image
set(gca, 'fontsize', 14)
xlim([0 1])
ylim([0 0.8])
xticks(linspace(0,1,11))
box on 

min_white_margin(gca);

%saveas(fig, 'point_sources', 'png')

merged = find(~cellfun(@isempty,selected));
area_vec = 1:length(merged);
for i = 1:length(merged)
    area_vec(i) = sum(cell_area(init_region{merged(i)}));
end

%(sum(area_vec<lbd)+1)/(length(area_vec)+1)
%figure
 %histogram(area_vec,100,'BinLimits',[0,0.02])