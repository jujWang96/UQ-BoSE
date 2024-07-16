load('real_full_result_2018_11_1_23_4_37.mat')
load('real_full_result_boot.mat')

GRAY = [0.6 0.6 0.6];

n_boot = size(impute_log_int, 2);

boot_std = std(impute_log_int, 0, 2);

fig = figure;
triplot(DT, 'Color', GRAY)
hold on
scatter(cx, cy, 12, log(boot_std), 'filled')
colorbar
colormap(jet)
axis image
min_white_margin(gca);
saveas(fig, 'log_boot_std', 'png')

fig = figure;
unique_n_region = unique(n_region);
counts = histc(n_region, unique_n_region);
% use percentage instead of count
bar(unique_n_region, counts / sum(counts), 'FaceColor', GRAY)
xlabel('Number of segments')
ylabel('Percentage')
set(gca, 'FontSize', 16)
min_white_margin(gca, 0, 0, 0.01);
