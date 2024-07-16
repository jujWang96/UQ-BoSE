function minindex = findminBIC(X,f_coef,M)
%find the number of fourier coefficients corresponds to minimum AIC
%Input:
% M: The maximum number of coefficients considered 

if nargin == 2
    M = 10;
end
Fbic = zeros(M,1);
for m = 1:M
    Fbic(m) = FourierBIC(X,f_coef,m);
end
[minBIC,minindex] = min(Fbic);
end