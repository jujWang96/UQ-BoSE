function  bic = FourierBIC(X,coef, M)
%
%X: input dataset
%coef: Fourier coefficients of segmentation contour
%M: number of Fourier coefficients used to recover the shape
%BIC = k*ln(n)-2*ln(L)
[nx,~] = size(X); 
[xv,yv] = iFD(coef,M);
%rmpath('folderthatisnotonpath')
% 
pgon = polyshape(xv,yv);
% 
if numboundaries(pgon)>1
    bic=0;
    return 
end
%            

A = polyarea(xv,yv);
in = inpolygon(X(:,1),X(:,2),xv,yv);
out = ~in;
bic = -2*(sum(in)*log(sum(in)/A)+sum(out)*log(sum(out)/(1-A)))+ 2*log(nx)*M;





