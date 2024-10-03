function CR = findCR(bootstrap_coef,num_coef,r,fileout)
%create the credible region defined in https://arxiv.org/abs/2007.04386
%input:
%   bootstap_coef,num_coef: used to parametries the the boundaris of
%   bootstrapped data, bounded in [0 1]x[0 1]
%   r: number of grid in row and column 
%output:
%   CI: a rxr matrix represent the binaries grids. 1 if the grid is in side
%   the 1-a credible region; 0 otherwise. 
%noteL If the pixel's central subpixel is inside the modified polygon, the pixel is inside the region.

hist = zeros(r,r);
n = length(num_coef);
for l = 1:n
    [invx,invy] = iFD(bootstrap_coef(:,l),num_coef(l));
    hist = hist+double(poly2mask(invx*r,invy*r,r,r));
end
%CR = zeros(r,r);
%CR(hist<(1-a/2)*n & hist>a/2*n)=1;

%CR = flip(CR);

CR = flip(hist/n);
CR = rot90(CR,3);
figure
histogram2('XBinEdges',0:r,'YBinEdges',0:r,'BinCounts',CR,'DisplayStyle','tile')
axis equal
colorbar
saveas(gcf,strcat(fileout,'.png'))


