function histogram_plot_group(n_region, true_n_source, factors, sample_factors)
% plot histograms of n_source

subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.06], [0.3 0.02], [0.08 0.02]);
% subtract the background region
%n_source = n_region - 1;
n_source = n_region;
index = 0;
for i = 1:length(factors)
    factor = factors(i);
    % only add xlabel to the last row
    if i == length(factors)
        add_xlabel = true;
    else
        add_xlabel = false;
    end
    for j = 1:length(sample_factors)
        index = index + 1;
        sample_factor = sample_factors(j);
        subplot_handle = subplot(length(factors), length(sample_factors), index);
        current_n_source = n_source(index, :); 
        histogram_plot(subplot_handle, current_n_source, true_n_source,...
            factor, sample_factor, add_xlabel)
        ax = gca;
        ax.FontSize = 15;
    end

end

end