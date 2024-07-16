function invcoef = make_scale_inv(coef)
%Normalize the FD so that |G|_2 = 1 
%input:
%  coef: original fourier coefficients
%output:
%  invcoef: fourier coefficients that is invariant to scale
    n = norm(coef(2));
    invcoef = coef/n;
end
