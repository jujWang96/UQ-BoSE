function [] = make_plot_GCR(X,origin_coef,boot_F_coef,gamma,min_modelselect,n_t,r)
rate = 2;
CI_combine = 0;
for i = 1:length(gamma)
    CI_combine = CI_combine+1-tubeCI(origin_coef,boot_F_coef,n_t,min_modelselect,chi2cdf(gamma(i)^2,2),r)/rate;
end
CI_combine = CI_combine/length(gamma);
rgb = cat(3,  ones(r,r), ones(r,r)*0.5,CI_combine*1.7);
for row = 1:r
    for col = 1:r
        if CI_combine(row,col)==1
            rgb(row,col,:) = 1;
        
        end
    end
end

%imshow(flip(rgb),'AlphaData', thisalpha)
imshow(flip(rgb))
[cx,cy] = iFD(origin_coef);
hold on
plot(cx*r,(1-cy)*r,'r','LineWidth',0.7);
scatter(X(:,1)*r,(1-X(:,2))*r,0.6,'K')
h = gca;
h.Visible = 'On';
set(gca,'XTick',[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]*r, 'YTick', [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]*r)
xticklabels({'0','0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1'})
yticklabels({'1','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'})
axis equal
axis([0 r 0 r])
box on
end
