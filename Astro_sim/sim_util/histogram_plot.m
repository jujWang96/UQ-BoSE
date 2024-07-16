function [] = histogram_plot(subplot_handle, current_n_source, true_n_source,...
    factor, sample_factor, add_xlabel)
% plot histogram of n_source

GRAY = [0.6 0.6 0.6];
% see http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html#2
RED = [0.8500, 0.3250, 0.0980];

% it is by default sorted
unique_n_source = unique(current_n_source);
counts = histc(current_n_source, unique_n_source);
% use percentage instead of count
bar(unique_n_source, counts / sum(counts), 'FaceColor', GRAY)
% highlight the bar with correct number
which_bar = find(unique_n_source == true_n_source, 1);
hold on
bar(unique_n_source(which_bar), counts(which_bar) / sum(counts), 'FaceColor', RED)
if add_xlabel
    xlabel('K')
end
ylim([0 1])

subplot_pos = get(subplot_handle, 'position');
add_annotation([0.7 0.7 0.08 0.08], subplot_pos, ['\beta=', num2str(sample_factor)])
add_annotation([0.7 0.5 0.08 0.08], subplot_pos, ['\sigma=', num2str(factor)])

end