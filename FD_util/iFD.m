function [x,y] = iFD(coef,M,L)
%
% Inverse the coefficient to actual boundary
%Input: 
% coef: all Fourier coefficients of a boundary
% M: the number of Fourier coefficients that will be used to reconstruct
% boundary
% L: the extended length
%Output:
% x,y: coordinates of reconstructed shape 
% L:
%
%Algorithm refer to principle of digital image processing algorithm 6.4 
%
n = length(coef);
if nargin == 1 
   M = n;
   L = n;

end
if nargin == 2 
   L = n;
end

M = min(M,n);
%x = zeros(N,1);
%y = zeros(N,1);


   
%for k = 1:N
 %   t = (k-1)/N;
  %  [x(k),y(k)] = singleiFD(coef,Mp,t);
    
%end



coef = [coef(1:ceil((M-1)/2)+1);zeros(L-M,1);coef((n-floor((M-1)/2)+1):n)];

a = ifft(coef*length(coef));

x = real(a);
y = imag(a);
output = [x,y];
end

