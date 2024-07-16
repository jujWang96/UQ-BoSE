function aicc = FourierAICc(X,coef,M)
%modify AIC for small sample size, AICc = AIC+(2k^2+2k)/(n-k-1)
[nx,~] = size(X); 
cx = X(:,1);
cy = X(:,2);
[xv,yv] = iFD(coef,M);
A = polyarea(xv,yv);
in = inpolygon(cx,cy,xv,yv);
out = ~in;
aicc = -2*(sum(in)*log(sum(in)/A)+sum(out)*log(sum(out)/(1-A)))+ 2*2*M+(2*(2*M)^2+2*2*M)/(nx-2*M-1);


