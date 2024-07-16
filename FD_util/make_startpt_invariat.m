function f_coef = make_startpt_inv(f_coef,Mp)
%make the contour sample start point invariate
%Mp: pair of cofficients used to make the sample invariate

N = numel(f_coef);
if nargin==1
    Mp = floor((N-1)/2);
end

fun = @(x) fphi(x,f_coef,Mp);


%find the phi that maximizes fphi in [0,pi) by brute-force search

cmax= -Inf;
phimax = 0;
A= 400;
for a = (-A/2):1:(A/2)
    phi = pi*(a-1)/400;
    c = fun(phi);
    if c>cmax
        cmax = c;
        phimax = phi;
    end
    
end

disp(phimax);

%shifting the starting point phase

for m =1:Mp
    f_coef(1+m) = exp(1i*m*phimax)*f_coef(1+m);
    f_coef(N-m+1) = exp(-1i*m*phimax)*f_coef(N-m+1);
end

%transform the entire coefficients
if Mp==floor((N-1)/2)
    if Mp<ceil((N-1)/2)
        f_coef(1+Mp+1) = exp(1i*(Mp+1)*phimax)*f_coef(1+Mp+1);
    end
else
end






