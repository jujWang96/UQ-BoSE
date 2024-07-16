function invcoef = make_rotate_inv(coef)
%input:
%  coef: original fourier coefficients
%output:
%  invcoef: fourier coefficients that is invariant to translation
    z = coef(2)+coef(end);
    b = atan(imag(z)/real(z));
    invcoef = coef;
    invcoef(2:end) = coef(2:end)*exp(-1i*b);
end
