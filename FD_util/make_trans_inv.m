function invcoef = make_trans_inv(coef)
%input:
%  coef: original fourier coefficients
%output:
%  invcoef: fourier coefficients that is invariant to translation
    invcoef = coef;
    invcoef(1) = 0;
end
