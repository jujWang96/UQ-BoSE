function histogram_plot_group_model_selection_singleshape(FD_nums, true_FD_num,factors, sample_factors)
% plot histograms of FD_num

subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.06], [0.1 0.02], [0.08 0.02]);

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
        current_FD_nums = FD_nums{index}; 
        histogram_plot_model_selection(subplot_handle, current_FD_nums, true_FD_num,...
            factor, sample_factor, add_xlabel,true)
        ax = gca;
        ax.FontSize = 15;
    end

end

end

