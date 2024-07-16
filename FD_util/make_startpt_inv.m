function [f_coef_a,f_coef_b] = make_startpt_inv(f_coef,Mp)
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
for a = 1:A
    phi = pi*(a-1)/400;
    c = fun(phi);
    if c>cmax
        cmax = c;
        phimax = phi;
    end
    
end


%shifting the starting point phase
f_coef_a = f_coef;
f_coef_b = f_coef;
for m =1:Mp
    f_coef_a(1+m) = exp(1i*m*phimax)*f_coef_a(1+m);
    f_coef_a(N-m+1) = exp(-1i*m*phimax)*f_coef_a(N-m+1);
end

%transform the entire coefficients or make the rest be 0
if Mp==floor((N-1)/2)
    if Mp<ceil((N-1)/2)
        f_coef_a(1+Mp+1) = exp(1i*(Mp+1)*phimax)*f_coef_a(1+Mp+1);
    end
else
    for m = (Mp+2):(N-Mp)
        f_coef_a(m)=0;
    end
end 
phimax = phimax+pi;    
for m =1:Mp
    f_coef_b(1+m) = exp(1i*m*phimax)*f_coef_b(1+m);
    f_coef_b(N-m+1) = exp(-1i*m*phimax)*f_coef_b(N-m+1);
end

%transform the entire coefficients or make the rest be 0
if Mp==floor((N-1)/2)
    if Mp<ceil((N-1)/2)
        f_coef_b(1+Mp+1) = exp(1i*(Mp+1)*phimax)*f_coef_b(1+Mp+1);
    end
else
    for m = (Mp+2):(N-Mp)
        f_coef_b(m)=0;
    end
end
end






