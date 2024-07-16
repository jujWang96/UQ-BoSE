function coef = FD(x,y)
%
% Calculate the fourier coefficients of a geometric boundary 
%Input: 
% x: nx1 x coordintes of a geometric boundary in counter-clockwise manner
% y: nx1 y coordintes of a geometric boundary in counter-clockwise manner
%    
    s = x + 1i*y;

    M = length(s);
    
    coef = fft(s)/M;
    
end
