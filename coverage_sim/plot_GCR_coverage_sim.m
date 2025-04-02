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

cd('~/UQ-BoSE/coverage_sim')
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
shapenames = ["ellipse","bridge","spiral"];
V(:,1) = V_spiral(:,1);
V(:,2) = V_spiral(:,2);
factor = 30;
sample_factor = 1;
lambda =1000;
base_num_in_circle = 10 ;
totalN = (factor*base_num_in_circle+lambda)*sample_factor;


for shapename = shapenames
    load(strcat(parpath,resultpath, "coverage_illustration_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor), "_gridsize", num2str(gridsize), ...
            "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search),'.mat'))
    
    subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.01], [0.05 0.02]);

    figure
    [invx,invy] = iFD(coef,candidate_aic);
    subplot(1,3,1)
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

    subplot(1,3,2)
    plot(CI_g1_aic)
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

    subplot(1,3,3)
    plot(CI_g2_aic)
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
    set(gcf,'Position',[100 75 800 300])

    saveas(gcf,strcat(parpath,plotpath,   "coverage_GCRAIC_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor), "_gridsize", num2str(gridsize), ...
            "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search),'_new'),'epsc')
    subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.03], [0.05 -0.02], [0.05 0.02]);

    figure
    [invx,invy] = iFD(coef,candidate_bic);
    subplot(1,3,1)
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

    subplot(1,3,2)
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

    subplot(1,3,3)
    plot(CI_g2_bic)
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
    set(gcf,'Position',[100 75 800 300])

    saveas(gcf,strcat(parpath,plotpath,   "coverage_GCRBIC_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor), "_gridsize", num2str(gridsize), ...
            "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search),'_new'),'epsc')

%%
end
