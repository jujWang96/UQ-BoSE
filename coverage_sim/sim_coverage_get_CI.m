function [BM_AICs,boundaries_AICs,BM_BICs, boundaries_BICs,origin_coef_inv,min_num] = sim_coverage_get_CI(V,totalN, objN, seed,rep_itr,gridsize,seedsize,B, p)
%when p is a vector, BM and boundaries are cells corresponding to each p
%value
X = sim_inhomo_const_general([0 1], [0 1], totalN-objN,{V(:,1)},{V(:,2)}, objN,seed,false);
[contourx,contoury, origin_coef] = sim_fit_coverage(X, seed, true,gridsize,seedsize,rep_itr);
%check if the original segmentation is valid 
K= 300;
BM_AICs = cell(1,length(p));
boundaries_AICs = cell(1,length(p));
BM_BICs = cell(1,length(p));
boundaries_BICs = cell(1,length(p));
origin_coef_inv = make_startpt_inv2(origin_coef);
min_num = cell(1,2);
if sum(contourx)==0
    return
end
out_bd_x = [0,1,1,0]';
out_bd_y = [0,0,1,1]';
bd = [41];
target_x{1} = contourx;
target_y{1} = contoury;
fourier_coef_all_AIC = zeros(K,B);
fourier_coef_all_BIC = zeros(K,B);

%[smooth_polyset,min_pair_num] = model_selection_multiple_fast(X,target_x,target_y,[{1},{1}],out_bd_x,out_bd_y,["AIC","BIC"],bd);

[smooth_polyset,min_num] = model_selection_multiple_fast_single(X,target_x,target_y,[{1},{1}],out_bd_x,out_bd_y,["AIC","BIC"],bd);
smooth_polyset_AIC = smooth_polyset{1};
parfor rep = 1:B
    bootstrap_X_set_AIC = sim_fix_data_multiple(X,[0,1], [0,1],smooth_polyset_AIC,1,rep);
    [~,~,fourier_coef_all_AIC(:,rep) ] = sim_fit_coverage(bootstrap_X_set_AIC{1}, rep, false,gridsize,seedsize,rep_itr);
    
end

fourier_coef_all_AIC( :, ~any(fourier_coef_all_AIC,1) ) = []; 

for ii = 1:length(p)
    [BM_AICs{ii},~,boundaries_AICs{ii}] = tubeCI_single_boundary(origin_coef,fourier_coef_all_AIC,300,min_num{1},p(ii),1000);
end
smooth_polyset_BIC = smooth_polyset{2};

parfor rep = 1:B
    bootstrap_X_set_BIC = sim_fix_data_multiple(X,[0,1], [0,1],smooth_polyset_BIC,1,rep);
    [~,~, fourier_coef_all_BIC(:,rep)] = sim_fit_coverage(bootstrap_X_set_BIC{1}, rep, false,gridsize,seedsize,rep_itr);
    
end
fourier_coef_all_BIC( :, ~any(fourier_coef_all_BIC,1) ) = []; 


for ii = 1:length(p)

    [BM_BICs{ii},~,boundaries_BICs{ii}] = tubeCI_single_boundary(origin_coef,fourier_coef_all_BIC,300,min_num{2},p(ii),1000);

end



