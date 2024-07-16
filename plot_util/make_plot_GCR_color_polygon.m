function [] = make_plot_GCR_color_polygon(X,origin_coef,boot_F_coef,gamma,min_modelselect,n_t,colors)
    r= 1024; 
    if nargin ==7
        polymask = false;
        for i = length(gamma):-1:1
            [CI,~,~] = tubeCI_single_boundary(origin_coef,boot_F_coef,n_t,min_modelselect,chi2cdf(gamma(i)^2,2),r, polymask);
            plot(CI,'FaceColor',colors(i,:),'FaceAlpha',1)
            hold on
            %rgb_combine(row,col,:) = colors(i,:);
        end
    else
        polymask = false;
        for i = length(gamma):-1:1
            [CI,~,~] = tubeCI_single_boundary(origin_coef,boot_F_coef,n_t,min_modelselect,chi2cdf(gamma(i)^2,2),r, polymask);
            plot(CI,'FaceAlpha',1)
            hold on
            %rgb_combine(row,col,:) = colors(i,:);
        end
    end
    [cx,cy] = iFD(origin_coef);
    plot([cx;cx(1)],([cy;cy(1)]),'r','LineWidth',2);
    scatter(X(:,1),(X(:,2)),1,'k.')
    h = gca;
    h.Visible = 'On';
    box on 
    axis equal
    axis([0 1 0 1])
    set(gcf,'renderer','Painters')
end