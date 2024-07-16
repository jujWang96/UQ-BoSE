addpath(genpath('~/Desktop/gsrg/NGC2300XMM'))
addpath(genpath('~/Desktop/gsrg/gcr_util'))
addpath(genpath('~/Desktop/gsrg/FD_util'))
addpath(genpath('~/Desktop/gsrg/gcr_util'))
addpath(genpath('~/Desktop/gsrg/plot_util'))
addpath(genpath('~/Desktop/gsrg/contour_util'))

addpath(genpath('~/Desktop/gsrg/Astro_sim/util'))
addpath(genpath('~/Desktop/gsrg/G-SRG'))

 %create scatter plot
load ngc2300_MOS1_evt_0.5-8.0keV_scaled.mat
calc_ra = @(x) calc_sec(x,4800.0,23707.5,26107.5,0.05,5.38516,1);
calc_dec = @(x) calc_sec(x,4800.0,25646.5,28046.5,0.05,-6.92139,2);
fig = figure;
X = double(unique(X,'rows'));
scatter(X(:,1),X(:,2),1,'K','.')
axis equal
axis image
box on
set(gca,'XTick',[calc_ra(100),calc_ra(50),calc_ra(0),calc_ra(-50),calc_ra(-100)],...
    'YTick', [calc_dec(-100),calc_dec(-50),calc_dec(0),calc_dec(50),calc_dec(100)],'Fontsize',18,'FontWeight','bold')
xticklabels({'100','50','0','-50','-100'})
yticklabels({'-100','-50','0','50','100'})
%set(gca, 'XDir','reverse')
%axis([calc_ra(60) calc_ra(-60) calc_dec(-60) calc_dec(60)])
axis([0 1 0 1 ])
xlabel('\DeltaRA [arcsec] +7:32:20.8511')
ylabel('\DeltaDec [arcsec] +85:42:32.186')

set(gcf, 'Position', [50 50 400 400]); %
set(gca,'XTick',[calc_ra(50),calc_ra(40),calc_ra(30),calc_ra(20),calc_ra(10),...
    calc_ra(0),calc_ra(-10),calc_ra(-20),calc_ra(-30),calc_ra(-40),calc_ra(-50)],...
    'YTick', [calc_dec(-50),calc_dec(-40),calc_dec(-30),calc_dec(-20),calc_dec(-10),...
    calc_dec(0),calc_dec(10),calc_dec(20),calc_dec(30),calc_dec(40),calc_dec(50)])
xticklabels({'50','40','30','20','10','0','-10','-20','-30','-40','-50'})
yticklabels({'-50','-40','-30','-20','-10','0','10','20','30','40','50'})
xlim([calc_ra(60),calc_ra(-60)])
ylim([calc_dec(-60),calc_dec(60)])
text(0.55, 0.75, 'NGC2300/XMM','FontSize',18,'FontWeight','bold')
set(gca,'XTickLabelRotation', 0,'linewidth',2)
saveas(gca,'NGC2300XMM_scatterplot','png')


%plot seeds for p25s10
P = 25;
threshold = 10;
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
[V, R] = voronoiDiagram(DT);

% get seeds
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);
%[seeds, num] = get_seeds_sim_kde(cx, cy,invalid, adj_mat, threshold);
[seeds, num] = get_seeds_sim_voronoi(cx, cy,invalid, adj_mat,cell_area,P, threshold);
fig = figure;
plot_seeds(DT, cx, cy, seeds, [], [], lines(num), num, 0)
set(gca,'XTick',[calc_ra(100),calc_ra(50),calc_ra(0),calc_ra(-50),calc_ra(-100)],...
    'YTick', [calc_dec(-100),calc_dec(-50),calc_dec(0),calc_dec(50),calc_dec(100)],'Fontsize',18,'FontWeight','bold')
xticklabels({'100','50','0','-50','-100'})
yticklabels({'-100','-50','0','50','100'})
%set(gca, 'XDir','reverse')
%axis([calc_ra(60) calc_ra(-60) calc_dec(-60) calc_dec(60)])
axis([0 1 0 1 ])
xlabel('\DeltaRA [arcsec] +7:32:20.8511')
ylabel('\DeltaDec [arcsec] +85:42:32.186')

set(gcf, 'Position', [50 50 400 400]); %
set(gca,'XTick',[calc_ra(50),calc_ra(40),calc_ra(30),calc_ra(20),calc_ra(10),...
    calc_ra(0),calc_ra(-10),calc_ra(-20),calc_ra(-30),calc_ra(-40),calc_ra(-50)],...
    'YTick', [calc_dec(-50),calc_dec(-40),calc_dec(-30),calc_dec(-20),calc_dec(-10),...
    calc_dec(0),calc_dec(10),calc_dec(20),calc_dec(30),calc_dec(40),calc_dec(50)])
xticklabels({'50','40','30','20','10','0','-10','-20','-30','-40','-50'})
yticklabels({'-50','-40','-30','-20','-10','0','10','20','30','40','50'})
xlim([calc_ra(60),calc_ra(-60)])
ylim([calc_dec(-60),calc_dec(60)])
text(0.55, 0.75, 'NGC2300/XMM','FontSize',18,'FontWeight','bold')
set(gca,'XTickLabelRotation', 0,'linewidth',2)

saveas(gca,'NGC2300XMM_seeds','png')

%plot segmentation result
fig = figure;
load ngc2300MOS_60_voronoi_p25s10_randstep60fraction1.mat   
colors = lines(length(glb_contour.regionId));
colors([3,4],:) = colors([4,3],:);
plot_segmentation_wo_voronoi(DT, glb_selected, cx, cy,colors ,8,true)
set(gca,'XTick',[calc_ra(100),calc_ra(50),calc_ra(0),calc_ra(-50),calc_ra(-100)],...
    'YTick', [calc_dec(-100),calc_dec(-50),calc_dec(0),calc_dec(50),calc_dec(100)],'Fontsize',18,'FontWeight','bold')

xticklabels({'100','50','0','-50','-100'})
yticklabels({'-100','-50','0','50','100'})
%set(gca, 'XDir','reverse')
%axis([calc_ra(60) calc_ra(-60) calc_dec(-60) calc_dec(60)])
axis([0 1 0 1 ])
xlabel('\DeltaRA [arcsec] +7:32:20.8511')
ylabel('\DeltaDec [arcsec] +85:42:32.186')

set(gcf, 'Position', [50 50 400 400]); %
set(gca,'XTick',[calc_ra(50),calc_ra(40),calc_ra(30),calc_ra(20),calc_ra(10),...
    calc_ra(0),calc_ra(-10),calc_ra(-20),calc_ra(-30),calc_ra(-40),calc_ra(-50)],...
    'YTick', [calc_dec(-50),calc_dec(-40),calc_dec(-30),calc_dec(-20),calc_dec(-10),...
    calc_dec(0),calc_dec(10),calc_dec(20),calc_dec(30),calc_dec(40),calc_dec(50)])
xticklabels({'50','40','30','20','10','0','-10','-20','-30','-40','-50'})
yticklabels({'-50','-40','-30','-20','-10','0','10','20','30','40','50'})
xlim([calc_ra(60),calc_ra(-60)])
ylim([calc_dec(-60),calc_dec(60)])
text(0.55, 0.75, 'NGC2300/XMM','FontSize',18,'FontWeight','bold')
set(gca,'XTickLabelRotation', 0,'linewidth',2)
box on
saveas(fig,'NGC2300XMM_segmentation','png')


%%%%%plot contours as well as the target contour
fig = figure;
target =  setdiff(setdiff( glb_contour.contourV{16}, glb_contour.contourV{20},'row'), glb_contour.contourV{24},'row');

make_plot_contour(X,glb_contour,1,2,target)

set(gca,'XTick',[calc_ra(100),calc_ra(50),calc_ra(0),calc_ra(-50),calc_ra(-100)],...
    'YTick', [calc_dec(-100),calc_dec(-50),calc_dec(0),calc_dec(50),calc_dec(100)],'Fontsize',18,'FontWeight','bold')
xticklabels({'100','50','0','-50','-100'})
yticklabels({'-100','-50','0','50','100'})
%set(gca, 'XDir','reverse')
%axis([calc_ra(60) calc_ra(-60) calc_dec(-60) calc_dec(60)])
axis([0 1 0 1 ])
xlabel('\DeltaRA [arcsec] +7:32:20.8511')
ylabel('\DeltaDec [arcsec] +85:42:32.186')

set(gcf, 'Position', [50 50 400 400]); %
set(gca,'XTick',[calc_ra(50),calc_ra(40),calc_ra(30),calc_ra(20),calc_ra(10),...
    calc_ra(0),calc_ra(-10),calc_ra(-20),calc_ra(-30),calc_ra(-40),calc_ra(-50)],...
    'YTick', [calc_dec(-50),calc_dec(-40),calc_dec(-30),calc_dec(-20),calc_dec(-10),...
    calc_dec(0),calc_dec(10),calc_dec(20),calc_dec(30),calc_dec(40),calc_dec(50)])
xticklabels({'50','40','30','20','10','0','-10','-20','-30','-40','-50'})
yticklabels({'-50','-40','-30','-20','-10','0','10','20','30','40','50'})
xlim([calc_ra(60),calc_ra(-60)])
ylim([calc_dec(-60),calc_dec(60)])
text(0.55, 0.75, 'NGC2300/XMM','FontSize',18,'FontWeight','bold')
set(gca,'XTickLabelRotation', 0,'linewidth',2)
box on
saveas(gca,'NGC2300XMM_contours','png')

%%%%%%%%
% plot the model selections
%%%%%%%
load ngc2300MOS_60_model_selection_voronoi_p25s10_1.mat
fig = figure;
make_plot_modelselect(smooth_polyset,1);
xlabel('\DeltaRA [arcsec] +7:32:20.8511')
ylabel('\DeltaDec [arcsec] +85:42:32.186')

set(gcf, 'Position', [50 50 400 400]); %
set(gca,'XTick',[calc_ra(50),calc_ra(40),calc_ra(30),calc_ra(20),calc_ra(10),...
    calc_ra(0),calc_ra(-10),calc_ra(-20),calc_ra(-30),calc_ra(-40),calc_ra(-50)],...
    'YTick', [calc_dec(-50),calc_dec(-40),calc_dec(-30),calc_dec(-20),calc_dec(-10),...
    calc_dec(0),calc_dec(10),calc_dec(20),calc_dec(30),calc_dec(40),calc_dec(50)],'Fontsize',18,'FontWeight','bold')
xticklabels({'50','40','30','20','10','0','-10','-20','-30','-40','-50'})
yticklabels({'-50','-40','-30','-20','-10','0','10','20','30','40','50'})
xlim([calc_ra(60),calc_ra(-60)])
ylim([calc_dec(-60),calc_dec(60)])
text(0.55, 0.75, 'NGC2300/XMM','FontSize',18,'FontWeight','bold')
set(gca,'XTickLabelRotation', 0,'linewidth',2)
box on
saveas(gca, 'NGC2300XMM_model_selection_aic', 'png')


fig = figure;

make_plot_modelselect(smooth_polyset,2);    
xlabel('\DeltaRA [arcsec] +7:32:20.8511')
ylabel('\DeltaDec [arcsec] +85:42:32.186')

set(gcf, 'Position', [50 50 400 400]); %
set(gca,'XTick',[calc_ra(50),calc_ra(40),calc_ra(30),calc_ra(20),calc_ra(10),...
    calc_ra(0),calc_ra(-10),calc_ra(-20),calc_ra(-30),calc_ra(-40),calc_ra(-50)],...
    'YTick', [calc_dec(-50),calc_dec(-40),calc_dec(-30),calc_dec(-20),calc_dec(-10),...
    calc_dec(0),calc_dec(10),calc_dec(20),calc_dec(30),calc_dec(40),calc_dec(50)],'Fontsize',18,'FontWeight','bold')
xticklabels({'50','40','30','20','10','0','-10','-20','-30','-40','-50'})
yticklabels({'-50','-40','-30','-20','-10','0','10','20','30','40','50'})
xlim([calc_ra(60),calc_ra(-60)])
ylim([calc_dec(-60),calc_dec(60)])
text(0.55, 0.75, 'NGC2300/XMM','FontSize',18,'FontWeight','bold')
set(gca,'XTickLabelRotation', 0,'linewidth',2)
box on
saveas(gca, 'NGC2300XMM_model_selection_bic', 'png')

%%%%plot GCR (BIC)
fig = figure;
load ngc2300MOS_60_voronoi_p25s10_randstep60fraction1.mat   
load ngc2300MOS_60_model_selection_voronoi_p25s10_1.mat
index=0;
boot_F_coefs = [];
for seed_num = 1:10
    load(strcat('ngc2300MOS_60_voronoip25s10_bootstrap_numerical_analysis_bic_seed',num2str(seed_num),'size20_fraction1.mat'))
     boot_F_coefs = [boot_F_coefs,boot_F_coef]; 
end
 
  load ngc2300_MOS1_evt_0.5-8.0keV_scaled.mat

objId = setdiff(seg_contour.regionId,seg_contour.background);
clear target_curve
[target_curve(:,1),target_curve(:,2)] = get_curve(seg_contour.contourV{objId},false);    
[centx,centy] = centroid(polyshape(target_curve(:,1),target_curve(:,2)));
[x_sample,y_sample] = sample_curve(target_curve(:,1),target_curve(:,2),300,centx, false);
origin_coef = FD(x_sample,y_sample);
r = 1024;
make_plot_GCR_color(X,origin_coef,boot_F_coefs,[1,2],min_pair_num{2}(3),1000,r);
xlabel('\DeltaRA [arcsec] +7:32:20.8511')
ylabel('\DeltaDec [arcsec] +85:42:32.186')

set(gcf, 'Position', [50 50 400 400]); %
set(gca,'XTick',[calc_ra(50),calc_ra(40),calc_ra(30),calc_ra(20),calc_ra(10),...
    calc_ra(0),calc_ra(-10),calc_ra(-20),calc_ra(-30),calc_ra(-40),calc_ra(-50)]*r,...
    'YTick', flip(r-[calc_dec(-50),calc_dec(-40),calc_dec(-30),calc_dec(-20),calc_dec(-10),...
    calc_dec(0),calc_dec(10),calc_dec(20),calc_dec(30),calc_dec(40),calc_dec(50)]*r),'Fontsize',18,'FontWeight','bold')
set(gca,'TickDir','in')

xticklabels({'50','40','30','20','10','0','-10','-20','-30','-40','-50'})
yticklabels(flip({'-50','-40','-30','-20','-10','0','10','20','30','40','50'}))
xlim([calc_ra(60),calc_ra(-60)]*r)
ylim([1-calc_dec(60),1-calc_dec(-60)]*r)
text(0.55*r, (1-0.75)*r, 'NGC2300/XMM','FontSize',18,'FontWeight','bold')
box on
a = get(gca,'XTickLabel');  
set(gca,'XTickLabelRotation', 0,'linewidth',2)
set(gcf,'renderer','Painters')
saveas(gca, 'NGC2300XMM_GCR_bic', 'png')

%%%%plot GCR (AIC)
fig = figure;
load ngc2300MOS_60_voronoi_p25s10_randstep60fraction1.mat   
load ngc2300MOS_60_model_selection_voronoi_p25s10_1.mat
     index=0;
     boot_F_coefs = [];
     for seed_num = 1:10
        load(strcat('ngc2300MOS_60_voronoip25s10_bootstrap_numerical_analysis_aic_seed',num2str(seed_num),'size20_fraction1.mat'))
         boot_F_coefs = [boot_F_coefs,boot_F_coef]; 
     end
 
  load ngc2300_MOS1_evt_0.5-8.0keV_scaled.mat

objId = setdiff(seg_contour.regionId,seg_contour.background);
clear target_curve
[target_curve(:,1),target_curve(:,2)] = get_curve(seg_contour.contourV{objId},false);    
[centx,centy] = centroid(polyshape(target_curve(:,1),target_curve(:,2)));
[x_sample,y_sample] = sample_curve(target_curve(:,1),target_curve(:,2),300,centx, false);
origin_coef = FD(x_sample,y_sample);
r = 1024;
make_plot_GCR_color(X,origin_coef,boot_F_coefs,[1,2],min_pair_num{1}(3),1000,r);
xlabel('\DeltaRA [arcsec] +7:32:20.8511')
ylabel('\DeltaDec [arcsec] +85:42:32.186')

set(gcf, 'Position', [50 50 400 400]); %
set(gca,'XTick',[calc_ra(50),calc_ra(40),calc_ra(30),calc_ra(20),calc_ra(10),...
    calc_ra(0),calc_ra(-10),calc_ra(-20),calc_ra(-30),calc_ra(-40),calc_ra(-50)]*r,...
    'YTick', flip(r-[calc_dec(-50),calc_dec(-40),calc_dec(-30),calc_dec(-20),calc_dec(-10),...
    calc_dec(0),calc_dec(10),calc_dec(20),calc_dec(30),calc_dec(40),calc_dec(50)]*r),'Fontsize',18,'FontWeight','bold')
set(gca,'TickDir','in')

xticklabels({'50','40','30','20','10','0','-10','-20','-30','-40','-50'})
yticklabels(flip({'-50','-40','-30','-20','-10','0','10','20','30','40','50'}))
xlim([calc_ra(60),calc_ra(-60)]*r)
ylim([1-calc_dec(60),1-calc_dec(-60)]*r)
text(0.55*r, (1-0.75)*r, 'NGC2300/XMM','FontSize',18,'FontWeight','bold')
box on
a = get(gca,'XTickLabel');  
set(gca,'XTickLabelRotation', 0,'linewidth',2)
set(gcf,'renderer','Painters')
saveas(gca, 'NGC2300XMM_GCR_aic', 'png')

%plot bootstrap flux 

num_in_O_seq_fds_aic = [];
 num_in_O_seqs_aic = [];

 for seed_num = 1:10
    load(strcat('ngc2300MOS_60_voronoip25s10_bootstrap_numerical_analysis_aic_seed',num2str(seed_num),'size20_fraction1.mat'))
     num_in_O_seq_fds_aic = [num_in_O_seq_fds_aic,num_in_O_seq_fd]; 
     num_in_O_seqs_aic = [num_in_O_seqs_aic,num_in_O_seq]; 

 end

fig = figure;
flux_aic{1} = num_in_O_seqs_aic;
flux_aic{2} = num_in_O_seq_fds_aic;
subplot(2,1,1)
make_hist_stack(flux_aic,[num_in,num_in_fd],linspace(2500,4000,51),{'bootstrap flux (non-truncated)','bootstrap flux (AIC)','obeserved flux','observerd aic flux'},'northeast')
ylim([0,0.15])


 num_in_O_seq_fds_bic = [];
  num_in_O_seqs_bic = [];

 for seed_num = 1:10
    load(strcat('ngc2300MOS_60_voronoip25s10_bootstrap_numerical_analysis_bic_seed',num2str(seed_num),'size20_fraction1.mat'))
     num_in_O_seq_fds_bic = [num_in_O_seq_fds_bic,num_in_O_seq_fd]; 
     num_in_O_seqs_bic = [num_in_O_seqs_bic,num_in_O_seq]; 

 end

subplot(2,1,2)
flux_bic{1} = num_in_O_seqs_bic;
flux_bic{2} = num_in_O_seq_fds_bic;
make_hist_stack(flux_bic,[num_in,num_in_fd],linspace(2500,4000,51),{'bootstrap flux (non-truncated)','bootstrap flux (BIC)','obeserved flux','observerd bic flux'},'northeast');
ylim([0,0.15])
saveas(gca, 'NGC2300XMM_bootstrap_flux', 'png')


%plot area
fig = figure;

subplot(2,1,1)
area_seq_fds_aic = [];
area_seqs_aic = [];

for seed_num = 1:10
    load(strcat('ngc2300MOS_60_voronoip25s10_bootstrap_numerical_analysis_aic_seed',num2str(seed_num),'size20_fraction1.mat'))
     area_seq_fds_aic = [area_seq_fds_aic,area_seq_fd]; 
     area_seqs_aic = [area_seqs_aic,area_seq]; 

end

area_aic{1} = area_seqs_aic;
area_aic{2} = area_seq_fds_aic;
make_hist_stack(area_aic,[area,area_fd],linspace(0.03,0.18,51),{'bootstrap area (non-truncated)','bootstrap area (AIC)','obeserved area','observerd aic area'},'northeast')
ylim([0,0.3])
subplot(2,1,2)
 area_seq_fds_bic = [];
  area_seqs_bic = [];

for seed_num = 1:10
    load(strcat('ngc2300MOS_60_voronoip25s10_bootstrap_numerical_analysis_bic_seed',num2str(seed_num),'size20_fraction1.mat'))
     area_seq_fds_bic = [area_seq_fds_bic,area_seq_fd]; 
     area_seqs_bic = [area_seqs_bic,area_seq]; 

end

area_bic{1} = area_seqs_bic;
area_bic{2} = area_seq_fds_bic;
make_hist_stack(area_bic,[area,area_fd],linspace(0.03,0.18,51),{'bootstrap area (non-truncated)','bootstrap area (BIC)','obeserved area','observerd bic area'},'northeast')
ylim([0,0.3])
saveas(gca, 'NGC2300XMM_bootstrap_area', 'png')

