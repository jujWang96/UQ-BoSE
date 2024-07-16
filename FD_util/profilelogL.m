function plofL = profilelogL(X,coef,M,lin,lout,n0 ) 

%Calculate the profile log likelihood of fitting X to coef under the poisson
%process assumption
%M: number of pairs choosen by model coef

%If no intensity is given, the intensity is esitmated based on X and the
%contour

cx = X(:,1);
cy = X(:,2);
[invx,invy] = iFD(coef,M);
in1 = inpolygon(cx,cy,invx,invy);
out1 = ~in1;
A = polyarea(invx,invy);
if nargin==3
    n0 = length(X);
    lin = sum(in1)/A;   
    lout = sum(out1)/(1-A);
else
    lin = lin/n0*length(X);
    lout = lout/n0*length(X);
    n0 = length(X);
end
%plofL = sum(in0)*log(sum(in0)/A)+sum(out0)*log(sum(out0)/(1-A));
plofL = sum(in1)*log(lin)+sum(out1)*log(lout)-n0;
end

