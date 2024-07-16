%load('real_full_all_seeds_result_2018_12_1_22_57_33.mat')
load('real_full_all_seeds_result_2022_11_20_13_50_21.mat')
addpath(genpath('~/Desktop/gsrg/Astro_sim/util'))
fig = figure;
triplot(DT, 'Color', GRAY)
hold on

index = 0;
selected_nonempty = {};
for i = 1:length(selected)
    if ~isempty(selected{i})
        index = index + 1;
        selected_nonempty{index} = selected{i};
    end
end

for i = 1:length(selected_nonempty)-1
    for j = i+1:length(selected_nonempty)
        % sorted by area from the largest to the smallest
        if sum(cell_area(selected_nonempty{i})) < sum(cell_area(selected_nonempty{j}))
            tmp = selected_nonempty{i};
            selected_nonempty{i} = selected_nonempty{j};
            selected_nonempty{j} = tmp;
        end
    end
end

for i = 1:length(selected_nonempty)
    log_int = log(sum(exp(cell_log_intensity(selected_nonempty{i})).*cell_area(selected_nonempty{i}))/sum(cell_area(selected_nonempty{i})));
    scatter(cx(selected_nonempty{i}), cy(selected_nonempty{i}), 12, log_int*ones(length(selected_nonempty{i}), 1), 'filled')
end
colorbar('EastOutside')
colormap(hsv)
caxis([10 18])

axis image
set(gca, 'fontsize', 12)
min_white_margin(gca);
saveas(fig, 'segment_result_points_all_seeds', 'png')
