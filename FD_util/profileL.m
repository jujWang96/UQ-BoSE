function plofL = profileL(X,coef,M) 

%Calculate the rpofile likelihood of fitting X to coef under the poisson
%process assumption
%M: number of pairs choosen by model coef
%should reject if LR is small 

if nargin==2
    M=length(coef);
end

cx = X(:,1);
cy = X(:,2);
[invx,invy] = iFD(coef,M);
A = polyarea(invx,invy);
in = inpolygon(cx,cy,invx,invy);
out = ~in;

plofL = sum(in)*log(sum(in)/A)+sum(out)*log(sum(out)/(1-A));

end

