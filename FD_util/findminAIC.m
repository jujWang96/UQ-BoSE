function minindex = findminAIC(X,f_coef,M)
%find the number of fourier coefficients corresponds to minimum AIC
%Input:
% M: The maximum number of coefficients considered 

if nargin == 2
    M = 20;
end
Faic = zeros(M,1);
for m = 1:M
    Faic(m) = FourierAIC(X,f_coef,m);
end
[minAIC,minindex] = min(Faic);
end