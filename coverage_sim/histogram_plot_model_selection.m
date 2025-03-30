function [] = histogram_plot_model_selection(subplot_handle, current_FD_num, true_FD_num,...
    factor, sample_factor, add_xlabel,isBIC)
% plot histogram of n_source

GRAY = [0.6 0.6 0.6];
% see http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html#2
RED = [0.8500, 0.3250, 0.0980];

% it is by default sorted
unique_FD_num = unique(current_FD_num);
counts = histc(current_FD_num, unique_FD_num);
% highlight the bar with correct number
which_bar = find(unique_FD_num == true_FD_num, 1);
% use percentage instead of count
bar(unique_FD_num, counts/ sum(counts), 0.8,'FaceColor', GRAY)

hold on
bar(unique_FD_num(which_bar), counts(which_bar) / sum(counts),1.6, 'FaceColor', RED)
if add_xlabel
    xlabel('J')
end
ylim([0 1])

subplot_pos = get(subplot_handle, 'position');
%add_annotation([0.7 0.7 0.08 0.08], subplot_pos, ['\beta=', num2str(sample_factor)])
%add_annotation([0.7 0.5 0.08 0.08], subplot_pos, ['\sigma=', num2str(factor)])
if ~isBIC
    selected_bars = [5,15,25,35,45,55];
    
    set(gca, 'xticklabels',({''})); 
    xticks(selected_bars)

    set(gca, 'xticklabels',(selected_bars)); 

end
end