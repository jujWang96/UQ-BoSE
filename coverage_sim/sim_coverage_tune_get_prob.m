function [GCRs,covers,origin_coef_inv] = sim_coverage_tune_get_prob(V,totalN, objN, seed,gridsize,seedsize, grid_only, merge_method,force_merge, rep_itr, pair_search, selection_method, boot_method,B, p, true_FD_num)
X = sim_inhomo_const_general([0 1], [0 1], totalN-objN,{V(:,1)},{V(:,2)}, objN,seed,false);
[contourx,contoury,coef,raw_contour] = sim_fit_coverage_tune(X, seed, false,gridsize,seedsize, grid_only, merge_method,force_merge, rep_itr);
%check if the original segmentation is valid 
K= 300;
GCRs = cell(1,length(p));
boundaries = cell(1,length(p));
covers = zeros(1,length(p));
origin_coef_inv = zeros(size(coef));
if all(contourx == 0)
    return
end
origin_coef_inv = make_startpt_inv2(coef);


if pair_search
    candidates = 3:2:K ;
    [candidate_aic,candidate_bic,pgon_confine] = model_selection_flexible_robust(X,contourx, contoury, raw_contour,candidates,false);
else
    candidates = 2:K ;
    [candidate_aic,candidate_bic,pgon_confine] = model_selection_flexible_robust(X,contourx, contoury, raw_contour,candidates,false);
end
fourier_coef_all = sim_coverage_tune_get_boot_coef(X, boot_method, coef, candidate_aic, candidate_bic,B,gridsize,seedsize, ...
    grid_only, merge_method,force_merge, rep_itr,true_FD_num);
use_polymask = false;
for ii = 1:length(p)
    if selection_method == "AIC"
        [GCRs{ii},~,boundaries{ii}] = tubeCI_single_boundary(coef,fourier_coef_all,K,candidate_aic,p(ii),1000,use_polymask);
    elseif selection_method == "BIC"
        [GCRs{ii},~,boundaries{ii}] = tubeCI_single_boundary(coef,fourier_coef_all,K,candidate_bic,p(ii),1000,use_polymask);
    else
        %use the true number of FD
        [GCRs{ii},~,boundaries{ii}] = tubeCI_single_boundary(coef,fourier_coef_all,K,true_FD_num,p(ii),1000,use_polymask);
    end
%     figure
%     plot(GCRs{ii})
%     hold on 
%     plot(V(:,1),V(:,2))
%     axis([0,1,0,1])
    point_inside = isinterior(GCRs{ii}, V);
    all_inside = all(point_inside);
    if all_inside
        covers(ii) = 1;
    end
end


end