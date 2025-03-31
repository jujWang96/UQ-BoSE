

function [metrics, n_region,rad] = sim_bulleye(loc, radius, unit_a, SNR, lambdas, contrasts, seed,gridsize,seedsize,rep_itr)

metrics = zeros(length(lambdas) * length(contrasts), 1);
n_region = zeros(length(lambdas) * length(contrasts), 1);
rad = zeros(length(lambdas) * length(contrasts), 1);
index = 0;
for i = 1:length(lambdas)
    lambda = lambdas(i);
    for j = 1:length(contrasts)
        index = index + 1;
        contrast = contrasts(j);
    
        base_num_in_circle = calc_num_in_source([1,3],unit_a,contrast,lambda,SNR,2);
        
        % generate simulated data (inhomogeneous Poisson point process)
        X = sim_inhomo_const_ring([0 1], [0 1] , loc, radius,1, lambda, base_num_in_circle, seed*length(lambdas) * length(contrasts)+(i*3+j), false);
        %save('debug.mat')
        if lambda==100
            gridsize = 0.15;
            seedsize = 1;
        elseif lambda == 500
            gridsize = 0.1;
            seedsize = 3;
        else
            gridsize = 0.1;
            seedsize = 5;
        end
        [labeled_cells, pred_class_all, n_region(index),rad(index)] = sim_fit_random_bulleye(X, seed*length(lambdas) * length(contrasts)+(i*3+j),false, gridsize,seedsize,rep_itr);
        true_class_all = evaluate_points_bulleye(X, loc, radius);
        % only care about labeled cells
        true_class_all = true_class_all(labeled_cells);
        pred_class_all = pred_class_all(labeled_cells);

        % it is symmetric
        metrics(index) = rand_index(true_class_all, pred_class_all, 'adjusted');
    end
end

end