function f = fphi(x,f_coef,Mp)
N = numel(f_coef);
f = 0;
for m = 1:Mp
    a = exp(-1i*m*x)*f_coef(N-m+1);
    b = exp(1i*m*x)*f_coef(m+1);
    f = f+real(a)*imag(b)-real(b)*imag(a);
end
end
