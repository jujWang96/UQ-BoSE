function [] = plot_kde(X,dx,dy,level)
    X(:, 1) = X(:, 1) *dx;
    X(:,2) = X(:, 2) *dy;
    itv = [0,dx,0,dy];
    [pdf_values, pdf_x,bw] = ksdensity(X, 'Kernel', 'epanechnikov','Function','pdf');
    kde_joint = @(x, y) griddata(pdf_x(:,1), pdf_x(:,2),pdf_values, x, y);
   
    figure;
    zhandle = fcontour(kde_joint,itv,'Fill','on','LevelStep',level);
    %axis([0,1,0,1])
    %axis equal
