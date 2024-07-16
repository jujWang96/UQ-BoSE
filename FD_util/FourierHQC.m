function hqc = FourierHQC(X,coef,M)
%Find the HQC 
%X: input dataset
%coef: Fourier coefficients of segmentation contour
%M: number of Fourier coefficients used to recover the shape

[nx,~] = size(X); 
cx = X(:,1);
cy = X(:,2);
[xv,yv] = iFD(coef,M);
A = polyarea(xv,yv);
in = inpolygon(cx,cy,xv,yv);
out = ~in;
hqc = -2*(sum(in)*log(sum(in)/A)+sum(out)*log(sum(out)/(1-A)))+ 2*log(log(nx))*M;




