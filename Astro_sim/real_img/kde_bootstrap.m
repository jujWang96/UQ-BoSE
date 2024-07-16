close all
clear
load('model_select_result_2300chandra.mat')
load('ngc2300_box_058kev_evt2.mat')
addpath(genpath('~/Desktop/gsrg/boot_util'))
cd '/Users/joycewang/Desktop/gsrg/NGC2300'
target_idx = 4;
X = double(unique(X,'rows'));
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
gridx1 = 0:.005:1;
gridx2 = 0:.005:1;
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
rng('default')  % For reproducibility

figure
ksdensity(X,xi,"Bandwidth",bw)

figure
f = mvksdensity(X,X);
[pdf_values, pdf_x] = mvksdensity(X, xi, 'Kernel', 'epanechnikov');
pdf_values_grid = reshape(pdf_values,[length(gridx2),length(gridx1)])';
N=10000
vals=zeros(2,N);
for i=1:N
[vals(1,i),vals(2,i)]=pinky(gridx1,gridx2,pdf_values_grid);
end
figure
scatter(vals(1,:),vals(2,:),'.')

imagesc(gridx1,gridx2,pdf_values_grid);
Y = datasample(X, 500);
scaling_factors = arrayfun(@(x, y) kde(x, y) * rand(), Y(:, 1), Y(:, 2));
samples = Y .* scaling_factors;

[cdf_values, cdf_x] = mvksdensity(X, xi,'Kernel', 'epanechnikov','Function','cdf','support',[-0.001,-0.001;1.001,1.001]);
kde_joint = @(x, y) griddata(cdf_x(:,1), cdf_x(:,2),cdf_values, x, y);


%[cdf_values, cdf_x] = mvksdensity(X(:,1), 'Kernel', 'epanechnikov','Function','cdf');
kde_x = @(x) interp1(cdf_x((cdf_x(:,2)==1),1), cdf_values((cdf_x(:,2)==1)), x);

[pdf_values, pdf_x] = mvksdensity(X, xi,'Kernel', 'epanechnikov','Function','pdf','support',[-0.001,-0.001;1.001,1.001]);
kde_pdf_x = @(x) interp1(pdf_x((pdf_x(:,2)==1),1), pdf_values((pdf_x(:,2)==1)), x);
integral(@(x) kde_pdf_x(x),0,1 )
%kde_x = @(x) integral(@(y) kde_joint(x,y), 0,1,'ArrayValued',true);
%kde_yonx = @(x, y) griddata(cdf_x(:,1), cdf_x(:,2),cdf_values, x, y)./interp1(cdf_x(:,1), cdf_values, x);
kde_yonx = @(x, y)kde_joint(x,y)./kde_pdf_x(x);
Y = generate_samples_from_cdf(kde_x,kde_yonx,100);





Y = generate_samples_from_kde(X,0.005:.005:0.995,0.005:.005:0.995,n);


[x0,y0]=pinky(Xin,Yin,dist_in);

[pdf_values, pdf_x] = mvksdensity(X, 'Kernel', 'epanechnikov','Function','pdf');
kde_joint = @(x, y) griddata(pdf_x(:,1), pdf_x(:,2),pdf_values, x, y);
zhandle = fcontour(kde_joint,'Fill','on');
axis([0,1,0,1])
axis equal
gca = image_rescaling(gca,'NGC2300');
colorbar
generate_samples_mcmc(kde_joint, 1000, 500, 10)
