function samples = generate_samples_from_kde(X,cdf_values, cdf_x,gridx1,gridx2, n_samples,seed)

    % Input:
    % Px - Function handle for the marginal CDF Px(x)
    % Pxy - Function handle for the conditional CDF Pxy(y|x)
    % n_samples - Number of samples to generate

    % Output:
    % samples - Generated 2D samples
    rng(seed)
%     cdf_fjoint = @(x, y) griddata(cdf_x(:,1), cdf_x(:,2),cdf_values, x, y);
%     cdf_fx = @(x) interp1(cdf_x((cdf_x(:,2)==max(cdf_x(:,2))),1), cdf_values((cdf_x(:,2)==max(cdf_x(:,2)))), x);
     Fx_dist = makedist('PiecewiseLinear', 'x', [0,cdf_x((cdf_x(:,2)==max(cdf_x(:,2))),1)',1], 'Fx',  [0,cdf_values((cdf_x(:,2)==max(cdf_x(:,2))))',1]);
    rand_x = random(Fx_dist, 1, n_samples);
    rand_y = zeros(1,n_samples);
        
%   rand_x = sort(random(Fx_dist, 1, n_samples));
%     [x1,x2] = meshgrid(rand_x, gridx2);
%     x1 = x1(:);
%     x2 = x2(:);
%     xi = [x1 x2];
%     [ccdf_values, ccdf_x] = mvksdensity(X, xi,'Kernel', 'epanechnikov','Function','cdf','support',[-0.001,-0.001;1.001,1.001]);
% 
%     ccdf_values = ccdf_values/max(ccdf_values);
%     ccdf_values(1) = 0;

    for i = 1:n_samples
%         if i==1
%             cdf_y_values = ccdf_values(ccdf_x(:,1)==rand_x(i));
%         else
%             cdf_y_values = ccdf_values(ccdf_x(:,1)==rand_x(i))-ccdf_values(ccdf_x(:,1)==rand_x(i-1));
%         end
         [x1,x2] = meshgrid(rand_x(i), gridx2);
         x1 = x1(:);
         x2 = x2(:);
         xi = [x1 x2];
        [pdf_values, pdf_x] = mvksdensity(X, xi,'Kernel', 'epanechnikov','Function','pdf','support',[-0.001,-0.001;1.001,1.001]);
        cpdf_values = cumsum(pdf_values)*(gridx1(2)-gridx1(1));
        cpdf_values = cpdf_values/max(cpdf_values);
        temp = gridx2;
        [v, w] = unique( cpdf_values, 'stable' );
        cpdf_values = cpdf_values(w);

        temp = temp(w);
        temp(cpdf_values<10e-5)=[];

         cpdf_values(cpdf_values<10e-5) = [];
         temp(cpdf_values>1-10e-5)=[];

        cpdf_values(cpdf_values>1-10e-5) = [];
%         cdf_y_values = cdf_y_values/max(cdf_y_values);
%         temp = ccdf_x(ccdf_x(:,1)==rand_x(i),2);
%         [v, w] = unique( cdf_y_values, 'stable' );
%         cdf_y_values = cdf_y_values(w);
%         temp = temp(w);
%         temp(cdf_y_values<10e-5)=[];
% 
%         cdf_y_values(cdf_y_values<10e-5) = [];
%         temp(cdf_y_values>1-10e-5)=[];
% 
%         cdf_y_values(cdf_y_values>1-10e-5) = [];

        %cdf_given_x = pdf_values(ccdf_x(:,1)==rand_x(i));
        %Fygivenx_dist = makedist('PiecewiseLinear', 'x', [0,cdf_x(cdf_x(:,1)==rand_x(i),2)',1], 'Fx',  cdf_given_x);
        %Fygivenx_dist = makedist('PiecewiseLinear', 'x',[0, temp',1], 'Fx',  [0,cdf_y_values',1]);
        Fygivenx_dist = makedist('PiecewiseLinear', 'x',[0, temp,1], 'Fx',  [0,cpdf_values',1]);

        rand_y(i) = random(Fygivenx_dist, 1, 1);
    end
    samples = [rand_x',rand_y'];
    

end