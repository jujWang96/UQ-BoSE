function coef = make_startpt_inv2(coef)
%make the coefficient start point invariant so that the phase of G1=0
    theta = -atan(imag(coef(2))/real(coef(2)));
    if real(coef(2))<0
        theta=theta-pi;
    end
    for m =1:length(coef)       
        pair = [m-1,m-length(coef)-1];
        [Mp,idx] = min(abs(pair));
        coef(m) = coef(m)*exp(1i*theta*Mp*sign(pair(idx)));
    end
end
