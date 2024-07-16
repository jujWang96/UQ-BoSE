function boxplot_group(data, g1, g2, name)
% data, g1, g2 are column vectors
% name is the name of the plot

bh = boxplot(data, {g1, g2}, 'colorgroup', g2, 'factorgap', 5, 'factorseparator', 1,...
    'Labels', {'0.5', '1', '2', '0.5', '1', '2', '0.5', '1', '2'});
set(bh(6, :), 'LineWidth', 1.5)

set(gca, 'FontSize', 16)

annotation('textbox', [0.2 0.9 0.08 0.08], 'String', '\sigma=10', 'FitBoxToText', 'on',...
    'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
annotation('textbox', [0.475 0.9 0.08 0.08], 'String', '\sigma=20', 'FitBoxToText', 'on',...
    'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
annotation('textbox', [0.75 0.9 0.08 0.08], 'String', '\sigma=30', 'FitBoxToText', 'on',...
    'FontSize', 16, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
if ~isempty(name)
    title(name)
end

end