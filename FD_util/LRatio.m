function lbd = LRatio(X,coef0,coef1,M0,M1) 
%Calculate the likelihood ratio between models 0 and model 1 under the poisson
%process assumption
%LR = ln(L0/L1), should reject if LR is small 
%M0: number of pairs choosen by model0


if nargin==3
    M0 = length(coef0);
    M1 = length(coef1);
end

cx = X(:,1);
cy = X(:,2);
[xv0,yv0] = iFD(coef0,M0);
A0 = polyarea(xv0,yv0);
in0 = inpolygon(cx,cy,xv0,yv0);
out0 = ~in0;
[xv1,yv1] = iFD(coef1,M1);
A1 = polyarea(xv1,yv1);
in1 = inpolygon(cx,cy,xv1,yv1);
out1 = ~in1;
lbd = (sum(in0)*log(sum(in0)/A0)+sum(out0)*log(sum(out0)/(1-A0))-...
    (sum(in1)*log(sum(in1)/A1)+sum(out1)*log(sum(out1)/(1-A1))));

end

