function [metrics, n_region] = sim5(range_x, range_y, base_num_in_L_shape,...
    loc, radius, base_num_in_circle, factors, lambda, sample_factors, seed)

metrics = zeros(length(factors) * length(sample_factors), 1);
n_region = zeros(length(factors) * length(sample_factors), 1);
index = 0;
for i = 1:length(factors)
    factor = factors(i);
    for j = 1:length(sample_factors)
        index = index + 1;
        sample_factor = sample_factors(j);
        disp(['factor=', num2str(factor), ',sample_factor=', num2str(sample_factor), ',seed=', num2str(seed)])
        % generate simulated data (inhomogeneous Poisson point process)
        X = sim_inhomo_Pois_const_L_shape(range_x, range_y, sample_factor * factor * base_num_in_L_shape, seed);

        X = [X; sim_inhomo_Pois_const([0 1], [0 1], sample_factor * lambda,...
            loc, radius, sample_factor * factor * base_num_in_circle)];

        [labeled_cells, pred_class_all, n_region(index)] = sim_fit(X, factor, sample_factor, seed, false);

        true_class_all = evaluate_points_sim5(X, range_x, range_y, loc, radius);
        % only care about labeled cells
        true_class_all = true_class_all(labeled_cells);
        pred_class_all = pred_class_all(labeled_cells);

        % it is symmetric
        metrics(index) = rand_index(true_class_all, pred_class_all, 'adjusted');
    end
end

end