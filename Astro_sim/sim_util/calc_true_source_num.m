function calc_true_source_num(n_region, true_n_source, factors, sample_factors)
% plot histograms of n_source

% subtract the background region
%n_source = n_region - 1;
n_source = n_region;
index = 0;
for i = 1:length(factors)
    factor = factors(i);
  
    for j = 1:length(sample_factors)
        index = index + 1;
        sample_factor = sample_factors(j);
        current_n_source = n_source(index, :);

        true_per = length(find(current_n_source==true_n_source))/length(current_n_source);
        disp(strcat('for factor ', num2str(factor), ' and sample factor ',num2str(sample_factor),' correct segmentation: ',num2str(true_per)));
    end
end

end