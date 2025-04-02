clear 
close all
parpath = '~/UQ-BoSE/';
plotpath = 'coverage_sim/plots/';
resultpath = 'coverage_sim/results/';
addpath(genpath(strcat(parpath,'Astro_sim')))
addpath(genpath(strcat(parpath,'G-SRG')))
addpath(genpath(strcat(parpath,'plot_util')))
addpath(genpath(strcat(parpath,'contour_util')))
addpath(genpath(strcat(parpath,'FD_util')))
addpath(genpath(strcat(parpath,'boot_util')))
addpath(genpath(strcat(parpath,'gcr_util')))

cd('~/src/gsrg/coverage_sim')
load('shape_coordinates.mat')
linewidth = 1.5;
%% srgong parameter
gridsize = 0.1;
seedsize = 5;
grid_only = true;
force_merge = true;
pair_search = true;
show_plot = true;
merge_method = "random";
rep_itr = 1000;
M = 11;
r = 1000;

%% FD and bootstrap parameters
B = 200;
K = 300;

GRAY = [0.6 0.6 0.6];

Beta = 0.25:0.25:4;

p = 1-exp(-Beta.^2/2);

%% data generate parameter
shapename = 'ellipse';
V(:,1) = V_ellipse(:,1);
V(:,2) = V_ellipse(:,2);
factor = 30;
sample_factor = 1;
lambda =1000;
base_num_in_circle = 10 ;
totalN = (factor*base_num_in_circle+lambda)*sample_factor;

seed = 0;

rng(seed)

objN = binornd(totalN,sample_factor*factor * base_num_in_circle/totalN);
%objN = poissrnd(sample_factor*factor * base_num_in_circle);
backgroundN = totalN - objN; %lambda*alpha;
X = sim_inhomo_const_general([0 1], [0 1], backgroundN,{V(:,1)},{V(:,2)}, objN,seed,false);
[contourx,contoury,coef,raw_contour] = sim_fit_coverage_tune(X, seed, show_plot,gridsize,seedsize, grid_only, merge_method,force_merge, rep_itr);
set(gcf,'Position',[100 75 800 300])
subplot(1,3,1)
scatter(X(:,1),X(:,2),'k.')
hold on
plot(V(:,1),V(:,2),'Color',GRAY,'Linewidth',linewidth)
box on 
axis equal
axis([0,1,0,1])
set(get(gca,'ylabel'),'rotation',0)
xlabel(texlabel('X'))
ylabel(texlabel('Y'))
set(gca,'fontsize',14)
xticks(0:0.2:1)
yticks(0:0.2:1)
saveas(gcf,strcat(parpath,plotpath, "coverage_illustration_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor)  , "_gridsize", num2str(gridsize), ...
            "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search)),'epsc')

if pair_search
    candidates = 3:2:K ;
    [candidate_aic,candidate_bic,pgon_confine] = model_selection_flexible_robust(X,contourx, contoury, raw_contour,candidates,true);
else
    candidates = 3:K ;
    [candidate_aic,candidate_bic,pgon_confine] = model_selection_flexible_robust(X,contourx, contoury, raw_contour,candidates,true);
end
disp(strcat("aic model selection ", num2str(candidate_aic), " FDs"))
disp(strcat("bic model selection ", num2str(candidate_bic), " FDs"))
%compare the two model selection method
% out_bd_x = [0,1,1,0]';
% out_bd_y = [0,0,1,1]';
% bd = [41];
% target_x{1} = contourx;
% target_y{1} = contoury;
% 
% [smooth_polyset,min_num] = model_selection_multiple_fast_single(X,target_x,target_y,[{1},{1}],out_bd_x,out_bd_y,["AIC","BIC"],bd);
% 

figure
subplot(1,2,1)

[invx,invy] = iFD(coef,candidate_aic);
plot(invx,invy)

box on 
axis equal
axis([0,1,0,1])
xlabel(texlabel('X'))
ylabel(texlabel('Y'))
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',14)
xticks(0:0.2:1)
yticks(0:0.2:1)
 subplot(1,2,2)

[invx,invy] = iFD(coef,candidate_bic);
plot(invx,invy)

box on 
axis equal
axis([0,1,0,1])
xlabel(texlabel('X'))
ylabel(texlabel('Y'))
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',14)
xticks(0:0.2:1)
yticks(0:0.2:1)

% [contourx,contoury, origin_coef] = sim_fit_coverage_(X, seed, true,gridsize,seedsize,rep_itr, "greedy");
set(gcf,'Position',[100 75 800 300])
saveas(gcf,strcat(parpath,plotpath,"coverage_selection_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor)  , "_gridsize", num2str(gridsize), ...
            "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search)),'epsc')

%% generate B bootstrap data and derive their fourier cofficients 
fourier_coef_all_aic = sim_coverage_tune_get_boot_coef(X, "AIC", coef, candidate_aic, candidate_bic,B,gridsize,seedsize, ...
    grid_only, merge_method,force_merge, rep_itr);

[CI_g1_aic,C,boundaries_g1_aic] = tubeCI_single_boundary(coef,fourier_coef_all_aic,K,candidate_aic,p(4),r,false);
[CI_g2_aic,C,boundaries_g2_aic] = tubeCI_single_boundary(coef,fourier_coef_all_aic,K,candidate_aic,p(8),r,false);

fourier_coef_all_bic = sim_coverage_tune_get_boot_coef(X, "BIC", coef, candidate_aic, candidate_bic,B,gridsize,seedsize, ...
    grid_only, merge_method,force_merge, rep_itr);

[CI_g1_bic,C,boundaries_g1_bic] = tubeCI_single_boundary(coef,fourier_coef_all_bic,K,candidate_bic,p(4),r,false);
[CI_g2_bic,C,boundaries_g2_bic] = tubeCI_single_boundary(coef,fourier_coef_all_bic,K,candidate_bic,p(8),r,false);


%% plot the global confidence region for beta = 1, 2
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.05], [0.05 0.02]);

figure
[invx,invy] = iFD(coef,candidate_aic);
subplot(2,3,1)
scatter(X(:,1),X(:,2),'k.')
hold on
plot(invx,invy,'r','Linewidth',linewidth)
box on 
axis equal
axis([0,1,0,1])
xlabel(texlabel('X'))
ylabel(texlabel('Y'))
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',14)
xticks(0:0.2:1)
yticks(0:0.2:1)
subplot(2,3,2)
plot(CI_g1_aic)
hold on
plot([V(:,1);V(1,1)],[V(:,2);V(1,2)],'Color',GRAY,'Linewidth',linewidth)

box on 
axis equal
axis([0,1,0,1])
xlabel(texlabel('X'))
ylabel(texlabel('Y'))
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',14)
xticks(0:0.2:1)
yticks(0:0.2:1)
subplot(2,3,3)
plot(CI_g2_aic)
hold on
plot([V(:,1);V(1,1)],[V(:,2);V(1,2)],'Color',GRAY,'Linewidth',linewidth)

box on 
axis equal
axis([0,1,0,1])
xlabel(texlabel('X'))
ylabel(texlabel('Y'))
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',14)
xticks(0:0.2:1)
yticks(0:0.2:1)

[invx,invy] = iFD(coef,candidate_bic);

subplot(2,3,4)
scatter(X(:,1),X(:,2),'k.')
hold on
plot(invx,invy,'r','Linewidth',1.5)
box on 
axis equal
axis([0,1,0,1])
xlabel(texlabel('X'))
ylabel(texlabel('Y'))
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',14)
xticks(0:0.2:1)
yticks(0:0.2:1)
subplot(2,3,5)
plot(CI_g1_bic)
hold on
plot([V(:,1);V(1,1)],[V(:,2);V(1,2)],'Color',GRAY,'Linewidth',linewidth)

box on 
axis equal
axis([0,1,0,1])
xlabel(texlabel('X'))
ylabel(texlabel('Y'))
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',14)
xticks(0:0.2:1)
yticks(0:0.2:1)
subplot(2,3,6)
plot(CI_g2_bic)
hold on
plot([V(:,1);V(1,1)],[V(:,2);V(1,2)],'Color',GRAY,'Linewidth',2)

box on 
axis equal
axis([0,1,0,1])
xlabel(texlabel('X'))
ylabel(texlabel('Y'))
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',14)
xticks(0:0.2:1)
yticks(0:0.2:1)
set(gcf,'Position',[100 75 800 555])
saveas(gcf,strcat(parpath,plotpath,   "coverage_GCR_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor), "_gridsize", num2str(gridsize), ...
            "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search)),'epsc')

%%
figure
subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.03], [0.05 0.02], [0.05 0.02]);

subplot(1,4,1)
scatter(X(:,1),X(:,2),'k.')
box on 
axis equal
hold on
plot([V(:,1);V(1,1)],[V(:,2);V(1,2)],'Color',GRAY,'Linewidth',linewidth)
axis([0,1,0,1])
xlabel(texlabel('X'))
ylabel(texlabel('Y'))
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',14)
xticks(0:0.2:1)
yticks(0:0.2:1)

subplot(1,4,2)
plot(contourx,contoury,'r','Linewidth',linewidth)
box on 
axis equal
axis([0,1,0,1])
xlabel(texlabel('X'))
set(gca,'fontsize',14)
xticks(0:0.2:1)
yticks(0:0.2:1)

subplot(1,4,3)
[invx,invy] = iFD(coef,candidate_bic);
plot(invx,invy,'r','Linewidth',linewidth)
box on 
axis equal
axis([0,1,0,1])
xlabel(texlabel('X'))
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',14)
xticks(0:0.2:1)
yticks(0:0.2:1)

subplot(1,4,4)
plot(CI_g1_bic)
hold on
plot([V(:,1);V(1,1)],[V(:,2);V(1,2)],'Color',GRAY,'Linewidth',linewidth)
box on 
axis equal
axis([0,1,0,1])
xlabel(texlabel('X'))
set(get(gca,'ylabel'),'rotation',0)
set(gca,'fontsize',14)
xticks(0:0.2:1)
yticks(0:0.2:1)
set(gcf,'Position',[100 75 1100 300])
saveas(gcf,strcat(parpath,plotpath,   "coverage_combine_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor), "_gridsize", num2str(gridsize), ...
            "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search)),'epsc')
save(strcat(parpath,resultpath, "coverage_illustration_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor), "_gridsize", num2str(gridsize), ...
            "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search),'.mat'))

%%