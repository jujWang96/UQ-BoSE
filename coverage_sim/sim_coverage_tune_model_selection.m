function [candidate_aic,candidate_bic] = sim_coverage_tune_model_selection(V,totalN, objN, seed,gridsize,seedsize, grid_only, merge_method,force_merge, rep_itr, pair_search, robust)
%refactorize the sim_coverage_get_CI method
% selection_method: model selection method to dertermine the number of
% coefficients in 
% chang
%when p is a vector, BM and boundaries are cells corresponding to each p
%value
if nargin==11
    robust= true;
end
X = sim_inhomo_const_general([0 1], [0 1], totalN-objN,{V(:,1)},{V(:,2)}, objN,seed,false);
[contourx,contoury,coef,raw_contour] = sim_fit_coverage_tune(X, seed, false,gridsize,seedsize, grid_only, merge_method,force_merge, rep_itr);
%check if the original segmentation is valid 
K= 300;
candidate_aic = 0; candidate_bic = 0;
if all(contourx==0)
    return
end
if  robust
    if pair_search
        candidates = 3:2:K ;
        [candidate_aic,candidate_bic,pgon_confine] = model_selection_flexible_robust(X,contourx, contoury, raw_contour,candidates,false);
    else
        candidates = 3:K ;
        [candidate_aic,candidate_bic,pgon_confine] = model_selection_flexible_robust(X,contourx, contoury, raw_contour,candidates,false);
    end
    
else
    out_bd_x = [0,1,1,0]';
    out_bd_y = [0,0,1,1]';
    bd = [300];
    target_x{1} = contourx;
    target_y{1} = contoury;
    if pair_search
        
        [~,min_num] = model_selection_multiple_fast(X,target_x,target_y,[{1},{1}],out_bd_x,out_bd_y,["AIC","BIC"],bd);
        candidate_aic = min_num{1}*2+1;
        candidate_bic = min_num{2}*2+1;
    else
        [~,min_num] = model_selection_multiple_fast_single(X,target_x,target_y,[{1},{1}],out_bd_x,out_bd_y,["AIC","BIC"],bd);
        candidate_aic = min_num{1};
        candidate_bic = min_num{2};
    end
    
end
    [invx,invy] = iFD(coef,candidate_bic);
%     fig = figure;
%     plot(X(:,1),X(:,2),'k.')
%     hold on
%     plot(invx,invy)
%     plot(contourx,contoury)
%     saveas(fig,strcat('model_selection_debug_seed',num2str(seed),"_robust",num2str(robust)),"epsc")
end


