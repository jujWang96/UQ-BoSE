function [] = make_hist_with_range(data,m,edges,qrange)
%plot histogram of given and indicate m and standard deviation 
    histogram(data,'BinEdges',edges,'Normalization','probability')
    % Indicate those on the plot.
    xline(m,'Color','g','LineWidth',2);
    %xline(mu, 'Color', 'r', 'LineWidth', 2);
    xline(prctile(data,qrange(1)), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    xline(prctile(data,qrange(2)), 'Color', 'r', 'LineWidth', 2, 'LineStyle', '--');
    yl = ylim;
    sMean = sprintf(['observed |G1| = %.3f,\n (%d,%d)^{th} percentile = (%.3f,%.3f)'],...
        m, qrange(1),qrange(2),prctile(data,qrange(1)),prctile(data,qrange(2)));
% Position the text 90% of the way from bottom to top.
    text(0, 0.9*yl(2), sMean, 'Color', 'k', ...
	'FontWeight', 'bold', 'FontSize', 10, ...
	'EdgeColor', 'b');
end
