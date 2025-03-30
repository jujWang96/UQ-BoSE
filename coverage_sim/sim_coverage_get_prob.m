function [cover_AICs, cover_BICs,origin_coef_inv ] = sim_coverage_get_prob(V,totalN, objN, seed,rep_itr,gridsize,seedsize,B, p)
X = sim_inhomo_const_general([0 1], [0 1], totalN-objN,{V(:,1)},{V(:,2)}, objN,seed,false);

[BM_AICs,~,BM_BICs, ~,origin_coef_inv,~] = sim_coverage_get_CI(X, seed,rep_itr,gridsize,seedsize,B, p);

cover_AICs = zeros(1, length(p));
cover_BICs = zeros(1, length(p));
if isempty(BM_AICs{1}) 
    return
end
for ii = 1:length(p)
    if check_CI_cover(V,BM_AICs{ii})
         cover_AICs(ii) = 1;
    end
    if check_CI_cover(V,BM_BICs{ii})
         cover_BICs(ii) = 1;
    end
end