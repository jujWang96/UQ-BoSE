function x = dist_complex_stand(G1,G2,M,pair)
%measure the mutual standardized distance of two contours based on the first M of
%coefficients or some specfic combination of coef
%input: 
%   pair: logic parameter, if true then use pairs of coefficient specified by
%   M, other wise, use single coefficients specified by M
%Output:
%   x: the euclidean distance between two fourier spectrum

G1 = abs(G1);
G2 = abs(G2);

if nargin==2
    pair = false;
    M = 1:numel(G1);
    
end
if nargin==3
    pair = false;
    
end
if pair == true
    n = length(G1);
    G1(1) = 0;
    G2(1) = 0;

    G1 = [G1(1:M+1);G1(n-M+1:n)];
    G2 = [G2(1:M+1);G2(n-M+1:n)];
    %x = sum(abs(G1-G2).^2)^(1/2);
    diff = (G1-G2)./sqrt(G1)./sqrt(G2);
    x = norm(diff);
else
    diff = (G1-G2)./sqrt(G1)./sqrt(G2);
    x = norm(diff);
end
    

end
