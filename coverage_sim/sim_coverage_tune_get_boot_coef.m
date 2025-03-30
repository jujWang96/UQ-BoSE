function [fourier_coef_all] = sim_coverage_tune_get_boot_coef(X, boot_method, origin_coef, aic, bic,B,gridsize,seedsize, grid_only, merge_method,force_merge, rep_itr ,true_FD_num)
% get the bootstrap coefficients using appropriate bootstrap method 
% allow more ad hoc parameter to be tuned 
% tuning parameter added: 
%   boot_method: "AIC", "BIC", "nonparametric"
%   note: aic, bic may not be used if the boot_method is
%   nonparametric
%   grid_only: boolean, true: using only grid seeds, false: add local max 
%   merge_method: "random", "greedy", (the rep_itr param is not used under greedy merger)
%   force_merge: boolean, true: force the results to merge to two with beam
%   merger, false: return default and drop the simulation 
%   

%check if the original segmentation is valid 
K = length(origin_coef);
fourier_coef_all = zeros(K,B);
out_bd_x = [0,1,1,0]';
out_bd_y = [0,0,1,1]';
if boot_method == "AIC"
    [invx,invy] = iFD(origin_coef,aic);
    smooth_polyset = {polyshape(invx,invy),polyshape({invx;out_bd_x},{invy;out_bd_y})};
    parfor rep = 1:B
        bootstrap_X_set = sim_fix_data_multiple(X,[0,1], [0,1],smooth_polyset,1,rep);
        [~,~,fourier_coef_all(:,rep),~] = sim_fit_coverage_tune(bootstrap_X_set{1}, rep, false,gridsize,seedsize, ...
            grid_only, merge_method,force_merge, rep_itr);
    end
elseif boot_method == "BIC"
    [invx,invy] = iFD(origin_coef,bic);
    smooth_polyset = {polyshape(invx,invy),polyshape({invx;out_bd_x},{invy;out_bd_y})};
    parfor rep = 1:B
        bootstrap_X_set = sim_fix_data_multiple(X,[0,1], [0,1],smooth_polyset,1,rep);
        [~,~,fourier_coef_all(:,rep),~] = sim_fit_coverage_tune(bootstrap_X_set{1}, rep, false,gridsize,seedsize, ...
            grid_only, merge_method,force_merge, rep_itr);
    end
elseif boot_method == "true_FD_num"
     [invx,invy] = iFD(origin_coef,true_FD_num);
    smooth_polyset = {polyshape(invx,invy),polyshape({invx;out_bd_x},{invy;out_bd_y})};
    parfor rep = 1:B
        bootstrap_X_set = sim_fix_data_multiple(X,[0,1], [0,1],smooth_polyset,1,rep);
        [~,~,fourier_coef_all(:,rep),~] = sim_fit_coverage_tune(bootstrap_X_set{1}, rep, false,gridsize,seedsize, ...
            grid_only, merge_method,force_merge, rep_itr);
    end
elseif boot_method == "nonparametric"
    [cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
     parfor rep = 1:B
        bootstrap_X_set = sim_fix_data_vor_fast_vk(X,cell_area, n, DT,1,seed);
        [~,~,fourier_coef_all(:,rep),~] = sim_fit_coverage_tune(bootstrap_X_set{1}, rep, false,gridsize,seedsize, ...
            grid_only, merge_method,force_merge, rep_itr);
     end
else
    disp("specify a correct bootstrap method")
end
delete(gcp)
invalid = all(fourier_coef_all==0);
fourier_coef_all = fourier_coef_all( :, ~invalid );

