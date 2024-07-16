addpath(genpath('~/Desktop/gsrg/Astro_sim'))
addpath(genpath('~/Desktop/gsrg/plot_util'))

subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.04 0.1], [0.04 0.04]);

load X_circle.mat
subplot(3,3,1)
scatter(X(:,1),X(:,2),0.5,'.')
axis([0 1 0 1])

axis square;
%set(ax,'yticklabel',[])

subplot(3,3,2)
segment_graph_contour(X,true,false,'r',2,true,false,2);
title('gSRG results', 'FontSize', 15);
subplot(3,3,3)
segment_graph_contour(X,false,true,'r',2,true,false,2);
title('segment boundary', 'FontSize', 15);

load check_X_pois_ellipse_018.mat
subplot(3,3,4)
X_check = check_X_set{1};
scatter(X_check(:,1),X_check(:,2),0.5,'.')
axis equal
axis ([0 1 0 1])

subplot(3,3,5)
segment_graph_contour(X_check,true,false,'r',2,true,false,2);

subplot(3,3,6)
segment_graph_contour(X_check,false,true,'r',2,true,2);

load X_zshape.mat
subplot(3,3,7)
scatter(X(:,1),X(:,2),0.5,'.')
axis([0 1 0 1])
axis square;

subplot(3,3,8)
segment_graph_contour(X,true,false,'r',2,true,2);
subplot(3,3,9)
segment_graph_contour(X,false,true,'r',2,true,2);
set(gca,'XTick',[], 'YTick', [])
set(ax,'yticklabel',[])



subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);




load X_circle.mat
[x1,x2,y1,y2,glb_min_BIC] = gSRG_random(X,true,true, 'k',2,true,2,2,20,1000,6); 
h1 = subplot(3,4,1)
cla(h1)
subplot(3,4,1)
scatter(X(:,1),X(:,2),0.5,'.')
hold on
hold on 
theta = 0:0.01:2*pi;
xrad = 0.2;
yrad = 0.2;
x = xrad * cos(theta) + 0.5;
y = yrad*sin(theta)+0.5;
plot(x,y)
axis square
title('true boundary', 'FontSize', 15);
subplot(3,4,2)
[invx, invy] = iFD(origin_coef_rand);
plot([invx;invx(1)],[invy;invy(1)],'r')
axis([0 1 0 1])
axis square;
title('segment boundary', 'FontSize', 15);


subplot(3,4,3)
M = findminIC(X,origin_coef,'BIC',true);
[invx, invy] = iFD(origin_coef,M);

plot([invx;invx(1)],[invy;invy(1)],'r')
axis([0 1 0 1])
axis square;
title('BIC result','FontSize', 15)

subplot(3,4,4)
M = findminIC(X,origin_coef,'AIC',true)
[invx, invy] = iFD(origin_coef,M);
plot([invx;invx(1)],[invy;invy(1)],'r')
axis([0 1 0 1])
axis square;
title('AIC result','FontSize', 15)

h1 = subplot(3,4,5);

load check_X_pois_ellipse_018.mat
X = check_X_set{1};
origin_coef = check_F_coef(:,1);
subplot(3,4,5)
scatter(X(:,1),X(:,2),0.5,'.')
hold on 
theta = 0:0.01:2*pi;
xrad = 0.18;
yrad = 4/18;
x = xrad * cos(theta) + 0.5;
y = yrad*sin(theta)+0.5;
plot(x,y)
axis equal
axis square

subplot(3,4,6)
[invx, invy] = iFD(origin_coef);

plot([invx;invx(1)],[invy;invy(1)],'r')
axis([0 1 0 1])
axis square;

subplot(3,4,7)
M = findminIC(X,origin_coef,'BIC',true)
[invx, invy] = iFD(origin_coef,M);

plot([invx;invx(1)],[invy;invy(1)],'r')
axis([0 1 0 1])
axis square;

subplot(3,4,8)
M = findminIC(X,origin_coef,'AIC',true)
[invx, invy] = iFD(origin_coef,M);
plot([invx;invx(1)],[invy;invy(1)],'r')
axis([0 1 0 1])
axis square;
title('segmentation contour recovered by AICc(7 pair)')

load X_zshape.mat



subplot(3,4,9)
scatter(X(:,1),X(:,2),0.5,'.')
hold on 
axis equal
axis square
hold on
cx = [0.2 0.6 0.6 0.8 0.8 0.4 0.4 0.2 0.2]; 
cy = [0.2 0.2 0.6 0.6 0.8 0.8 0.4 0.4 0.2];
plot(cx,cy);

subplot(3,4,10)
[invx, invy] = iFD(origin_coef);

plot([invx;invx(1)],[invy;invy(1)],'r')
axis([0 1 0 1])
axis square;

subplot(3,4,11)
M = findminIC(X,origin_coef,'BIC',true,20)
[invx, invy] = iFD(origin_coef,M);

plot([invx;invx(1)],[invy;invy(1)],'r')
axis([0 1 0 1])
axis square;
title('segmentation contour recovered by BIC(4 pair)')

subplot(3,4,12)
M = findminIC(X,origin_coef,'AIC',true,40)
[invx, invy] = iFD(origin_coef,M);
plot([invx;invx(1)],[invy;invy(1)],'r')
axis([0 1 0 1])
axis square;
set(gca,'XTick',[], 'YTick', [])
set(ax,'yticklabel',[])

title('segmentation contour recovered by AICc(10 pair)')

subplot(2,3,1)
load X_circle.mat
scatter(X(:,1),X(:,2),12,'filled')

load check_X_pois_020.mat
subplot(2,3,2)
X_check = check_X_set{1};
scatter(X_check(:,1),X_check(:,2),12,'filled')

subplot(2,3,3)
X_check = check_X_set{2};
scatter(X_check(:,1),X_check(:,2),12,'filled')
subplot(3,4,4)
X_check = check_X_set{4};
scatter(X_check(:,1),X_check(:,2),12,'filled')

load boot_circle_1000_denoise.mat
subplot(2,3,4)
X_boot = bootstrap_X_set_denoise{1};
scatter(X_boot(:,1),X_boot(:,2),12,'filled')
title("bootstrap data1")
subplot(2,3,5)
X_boot = bootstrap_X_set_denoise{2};
scatter(X_boot(:,1),X_boot(:,2),12,'filled')
title("bootstrap data2")

subplot(2,3,6)
X_boot = bootstrap_X_set_denoise{3};
scatter(X_boot(:,1),X_boot(:,2),12,'filled')
title("bootstrap data3")


subplot(1,3,1)
load LRtest_radi.mat
plot((0.17:0.01:0.23),LRpower,'-o')
hold on
plot((0.17:0.01:0.23),LRpower_AIC,'-o')
xlabel("radius")
ylabel("rejection rate")
legend("no selection","AICc",'Location','southwest')
title("significant level=0.1, circle of different radius")
subplot(1,3,2)


load LRtest_center.mat
plot((0.5:0.01:0.54),LRpower,'-o')
hold on 
plot((0.5:0.01:0.54),LRpower_AIC,'-o')
xlabel("center")
ylabel("rejection rate")
legend("no selection","AICc",'Location','southeast')
title("significant level=0.1, circle of different center")

subplot(1,3,3)

load LRtest_ellipse.mat
plot((0.17:0.01:0.24),LRpower,'-o')

hold on
plot((0.17:0.01:0.24),LRpower_AIC,'-o')
xlabel("radius")
ylabel("rejection rate")
legend("no selection","AICc",'Location','southwest')
title("significant level=0.1, ellipse of different horizontal center")




subplot(3,4,8)
X_boot = bootstrap_X_set_denoise{4};
scatter(X_boot(:,1),X_boot(:,2),12,'filled')
subplot(3,4,9)
load check_X_pois_019.mat
X_check = check_X_set{1};
scatter(X_check(:,1),X_check(:,2),12,'filled')
title("circle with radius=0.19")
subplot(3,4,10)
load check_X_pois_021.mat
X_check = check_X_set{1};
scatter(X_check(:,1),X_check(:,2),12,'filled')
title("circle with radius=0.21")
subplot(3,4,11)
load check_X_pois_center_051.mat
X_check = check_X_set{1};
scatter(X_check(:,1),X_check(:,2),12,'filled')
title("circle with center = (0.51,0.5)")
subplot(3,4,12)
load check_X_pois_ellipse_021.mat
X_check = check_X_set{1};
scatter(X_check(:,1),X_check(:,2),12,'filled')
title("ellipse with horizontal radius=0.21")







load boot_circle_1000_denoise.mat
subplot(3,4,9)
X_boot = bootstrap_X_set_denoise{1};
scatter(X_boot(:,1),X_boot(:,2),12,'filled')

subplot(3,4,10)
X_boot = bootstrap_X_set_denoise{2};
scatter(X_boot(:,1),X_boot(:,2),12,'filled')
subplot(3,4,11)
X_boot = bootstrap_X_set_denoise{3};
scatter(X_boot(:,1),X_boot(:,2),12,'filled')
subplot(3,4,12)
X_boot = bootstrap_X_set_denoise{4};
scatter(X_boot(:,1),X_boot(:,2),12,'filled')


load /Users/joycewang/Desktop/segmentation/Astro_sim/virtual_result/bootstrap_500/boot_X_set_16001.mat
for j = 1:100
   X = check_X_set{j};
   scatter(X(:,1),X(:,2),0.1,'.','k')
   hold on 
end
axis square



theta =0:0.01:2*pi;
V(:,1) = 0.5+cos(theta)*0.2;
V(:,2) = 0.5+sin(theta)*0.2;

subplot(2,2,1)    
imshow(CI)
plot_circles([0.5 0.5]*r, 0.2*r,1.5);
subplot(2,2,2)
imshow(CI)
plot_circles([0.5 0.5]*r, 0.2*r,1.5);


subplot = @(m,n,p) subtightplot (m, n, p, [0.02 0.02], [0.10 0.06], [0.06 0.06]);

%idx = [1/1.4,1/1.3,1/1.2,1/1.1,1,1.1,1.2,1.3,1.4];
idx = -0.04:0.01:0.04;
%plot ct
subplot(1,3,1)
load ct.mat
var=1;
h = zeros(5,1);
h(1) = plot(idx,B3200(:,4),'k--'); hold on;
h(2) = plot(idx,B3200(:,4),'k-'); 
h(3) = plot(idx,B3200(:,4),'r');
h(4) = plot(idx,B3200(:,4),'b');
h(5) = plot(idx,B3200(:,4),'Color','#77AC30');
plot(idx,B3200(:,var),'r-.')
hold on
plot(idx,B3200(:,4),'r')
plot(idx,B2600(:,var),'b--')
plot(idx,B2600(:,4),'b')
plot(idx,B1600(:,var),'--','Color','#77AC30')
plot(idx,B1600(:,4),'Color','#77AC30')
ylim([0,1])
ylabel('Empirical rejection rate','Fontsize',12)
xticks(-0.04:0.01:0.04)
xticklabels({'-.04','-.03','-.02','-.01','0','.01','.02','.03','.04'})

xlabel(texlabel('c'),'Fontsize',12)
legend([h(1),h(2)],{texlabel('xi_C'),texlabel('xi_{LR}')},'Location','southwest')
title('(a)','Fontsize',15)

%%%%plot sz
subplot(1,3,2)
load sz.mat
var=2;
h = zeros(5,1);
h(1) = plot(idx,B3200(:,4),'k-.'); hold on;
h(2) = plot(idx,B3200(:,4),'k-'); 
h(3) = plot(idx,B3200(:,4),'r');
h(4) = plot(idx,B3200(:,4),'b');
h(5) = plot(idx,B3200(:,4),'g');
plot(idx,B3200(:,var),'r-.')
hold on
plot(idx,B3200(:,4),'r')
plot(idx,B2600(:,var),'b-.')
plot(idx,B2600(:,4),'b')
plot(idx,B1600(:,var),'g-.')
plot(idx,B1600(:,4),'g')
ylim([0,1])
xlabel(texlabel('d'),'Fontsize',12)
title('(b)','Fontsize',15)

legend([h(1),h(2)],{texlabel('xi_D'),texlabel('xi_{LR}')},'Location','southwest')
xticks(-0.04:0.01:0.04)
xticklabels({'-.04','-.03','-.02','-.01','0','.01','.02','.03','.04'})
ah1=axes('position',get(gca,'position'),'visible','off');
legend(ah1,[h(3),h(4),h(5)],{texlabel('nu=2,eta=60'),texlabel('nu=2,eta=30'),texlabel('nu=1,eta=60')},'Location','southeast')


%plot sp
load sp.mat
var=3;
subplot(1,3,3)
h = zeros(5,1);
h(1) = plot(idx,B1600(:,var),'k:','LineWidth',1.7); hold on;
h(2) = plot(idx,B3200(:,4),'k-'); 
h(3) = plot(idx,B3200(:,4),'r');
h(4) = plot(idx,B3200(:,4),'b');
h(5) = plot(idx,B3200(:,4),'g');
plot(idx,B3200(:,var),':','Color','r','LineWidth',1.7)
hold on
plot(idx,B3200(:,4),'r')
plot(idx,B2600(:,var),':','Color','b','LineWidth',1.7)
plot(idx,B2600(:,4),'b')
plot(idx,B1600(:,var),':','Color','g','LineWidth',1.7)
plot(idx,B1600(:,4),'g')
ylim([0,1])
xlabel(texlabel('r'),'Fontsize',12)
%ylabel('Empirical rejection rate','Fontsize',12)
%legend({texlabel('xi_{LR}'),texlabel('xi_c'),texlabel('xi_{LR}'),texlabel('xi_c'),texlabel('xi_{LR}'),texlabel('xi_c')},'Location','southwest')
xticks(-0.04:0.01:0.04)
xticklabels({'^{1}/_{1.4}','^{1}/_{1.3}','^{1}/_{1.2}','^{1}/_{1.1}','1','1.1','1.2','1.3','1.4'})
title('(c)','Fontsize',15)

%legend(h,{texlabel('xi_{S_{BIC}}'),texlabel('xi_{LR}'),texlabel('alpha=2,beta=60'),texlabel('alpha=2,beta=30'),texlabel('alpha=1,beta=60')},'Location','southeast')
legend([h(1),h(2)],{texlabel('xi_{S_{BIC}}'),texlabel('xi_{LR}')},'Location','southwest')
set(gcf, 'Position', [50 50 700 370]); %



subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.1 0.1], [0.06 0.06]);
CR15 = tubeCI(origin_coef,boot_F_coef,1000,1,chi2cdf(1.5^2,2),1000);
   CR25 = tubeCI(origin_coef,boot_F_coef,1000,1,chi2cdf(2.5^2,2),1000);
   
load X_zshape_boot.mat
CR15z = tubeCI(origin_coef,boot_F_coef,1000,10,chi2cdf(1.5^2,2),1000);
CR25z = tubeCI(origin_coef,boot_F_coef,1000,10,chi2cdf(2.5^2,2),1000);

subplot(2,2,1)
imshow(CR15)
r=1000;
    plot_circles([0.5 0.5]*r, 0.2*r,1);
    title(texlabel('gamma=1.5, M_p=1(BIC)'),'Fontsize',15)
 

subplot(2,2,2)
imshow(CR25)
    plot_circles([0.5 0.5]*r, 0.2*r,1);
 title(texlabel('gamma=2.5, M_p=1(BIC)'),'Fontsize',15)

subplot(2,2,3)
imshow(CR15z)
r=1000;
cx = [0.2 0.6 0.6 0.8 0.8 0.4 0.4 0.2 0.2]; 
cy = [0.2 0.2 0.6 0.6 0.8 0.8 0.4 0.4 0.2];
hold on
plot((1-cx)*r,cy*r,'r','LineWidth',1);
 title(texlabel('gamma=1.5, M_p=10(AIC)'),'Fontsize',15)

subplot(2,2,4)
imshow(CR25z)
cx = [0.2 0.6 0.6 0.8 0.8 0.4 0.4 0.2 0.2]; 
cy = [0.2 0.2 0.6 0.6 0.8 0.8 0.4 0.4 0.2];
hold on
plot((1-cx)*r,cy*r,'r','LineWidth',1);
 title(texlabel('gamma=2.5, M_p=10(AIC)'),'Fontsize',15)

 load acisf04614_repro_evt2_scaled.mat
scatter(X(:,1),X(:,2),0.5,'.');
axis square
 title('1987a-acisf04614','Fontsize',15)
set(gcf, 'Position', [50 50 500 500]); %

load ring_seed4_poisson.mat
load ring_seed4_greedy-025255.mat
fontsize=12;
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.02], [0.03 0.05], [0.03 0.03]);
subplot(3,3,1)
X = X_hh;
disp(length(X))
contour = contour_hh;
selected = selected_hh;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(double(X), [0 1], [0 1], ones(size(X, 1), 1));
plot_segmentation(DT, selected, cx, cy, lines(length(contour.regionId)))
box on
title(texlabel('(a)'),'Fontsize',fontsize)
subplot(3,3,2)
X = X_hm;
disp(length(X))
contour = contour_hm;
selected = selected_hm;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
plot_segmentation(DT, selected, cx, cy, lines(length(contour.regionId)))
box on
title(texlabel('(b)'),'Fontsize',fontsize)
subplot(3,3,3)
X = X_hl;
disp(length(X))

contour = contour_hl;
selected = selected_hl;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
plot_segmentation(DT, selected, cx, cy, lines(length(contour.regionId)))
box on
title(texlabel('(c)'),'Fontsize',fontsize)
subplot(3,3,4)
X = X_mh;
disp(length(X))
contour = contour_mh;
selected = selected_mh;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
plot_segmentation(DT, selected, cx, cy, lines(length(contour.regionId)))
box on
title(texlabel('(d)'),'Fontsize',fontsize)
subplot(3,3,5)
X = X_mm;
disp(length(X))
contour = contour_mm;
selected = selected_mm;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
plot_segmentation(DT, selected, cx, cy, lines(length(contour.regionId)))
box on
title(texlabel('(e)'),'Fontsize',fontsize)
subplot(3,3,6)
X = X_ml;
disp(length(X))
contour = contour_ml;
selected = selected_ml;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
plot_segmentation(DT, selected, cx, cy, lines(length(contour.regionId)))
box on
title(texlabel('(f)'),'Fontsize',fontsize)
subplot(3,3,7)
X = X_lh;
disp(length(X))
contour = contour_lh;
selected = selected_lh;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
plot_segmentation(DT, selected, cx, cy, lines(length(contour.regionId)))
box on
title(texlabel('(g)'),'Fontsize',fontsize)
subplot(3,3,8)
X = X_lm;
disp(length(X))
contour = contour_lm;
selected = selected_lm;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
plot_segmentation(DT, selected, cx, cy, lines(length(contour.regionId)))
box on
title(texlabel('(h)'),'Fontsize',fontsize)
subplot(3,3,9)
X = X_ll;
disp(length(X))
contour = contour_ll;
selected = selected_ll;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
plot_segmentation(DT, selected, cx, cy, lines(length(contour.regionId)))
box on
title(texlabel('(i)'),'Fontsize',fontsize)
set(gcf, 'Position', [50 50 600 600]); %


