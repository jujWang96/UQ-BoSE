function histogram_plot_group_model_selection(FD_nums, true_FD_nums,factor, sample_factor,isBIC)
% plot histograms of n_source

subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.06], [0.22 0.05], [0.08 0.02]);
% subtract the background region
%n_source = n_region - 1;
index = 0;
for i = 1:length(true_FD_nums)
    true_FD_num = true_FD_nums(i);
    add_xlabel = true;
   
    index = index + 1;
    subplot_handle = subplot( 1,length(true_FD_nums), index);
    curr_FD_num = FD_nums{index}; 
    histogram_plot_model_selection(subplot_handle, curr_FD_num, true_FD_num,...
        factor, sample_factor, add_xlabel,isBIC)
    ax = gca;
    ax.FontSize = 15;
    

end

end