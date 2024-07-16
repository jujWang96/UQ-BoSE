function [] = make_plot_GCR_color(X,origin_coef,boot_F_coef,gamma,min_modelselect,n_t,r,colors)
rgb_combine = cat(3,  ones(r,r), ones(r,r),ones(r,r));
if nargin == 7
    colors = [[0.8,1,0.8];[0.8,0.8,1]];
end
for i = length(gamma):-1:1
    CI = tubeCI(origin_coef,boot_F_coef,n_t,min_modelselect,chi2cdf(gamma(i)^2,2),r);
    [row,col ] = find(CI==1);
    for j = 1:length(row)
        
        rgb_combine(row(j),col(j),1) = colors(i,1);
        rgb_combine(row(j),col(j),2) = colors(i,2);

        rgb_combine(row(j),col(j),3) = colors(i,3);

        

    end
    
    %rgb_combine(row,col,:) = colors(i,:);
end



%imshow(flip(rgb),'AlphaData', thisalpha)
imshow(flip(rgb_combine),'InitialMagnification','fit')
[cx,cy] = iFD(origin_coef);
hold on
plot([cx;cx(1)]*r,(1-[cy;cy(1)])*r,'r','LineWidth',2);
scatter(X(:,1)*r,(1-X(:,2))*r,1,'k.')
h = gca;
h.Visible = 'On';

set(gcf,'renderer','Painters')
end
