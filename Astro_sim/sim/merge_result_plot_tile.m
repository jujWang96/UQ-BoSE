%load 
% 
close all;
addpath(genpath('~/src/gsrg/Astro_sim/sim'))
addpath(genpath('~/src/gsrg/plot_util'))

T  = 500;
gridsizes = [0.05,0.1,0.2];
factors = [10,20,30];
sample_factors = [0.5,1,2];
K = length(gridsizes);
metrics_all = zeros(K,length(factors)*length(sample_factors),T);
metrics_beam_all = zeros(K,length(factors)*length(sample_factors),T);
metrics_random_all = zeros(K,length(factors)*length(sample_factors),T);
metrics_all_fix = zeros(K,length(factors)*length(sample_factors),T);
metrics_beam_all_fix = zeros(K,length(factors)*length(sample_factors),T);
metrics_random_all_fix = zeros(K,length(factors)*length(sample_factors),T);

for k = flip(1:K)
    gridsize = gridsizes(k);
    load(strcat('sim4_result_gs',num2str(gridsize),'.mat'))
    load(strcat('sim4_result_random_gs',num2str(gridsize),'.mat'))
    load(strcat('sim4_result_beam_gs',num2str(gridsize),'.mat'))
%     if k==K
%         load(strcat('sim4_result_gridonly_gs',num2str(gridsize),'.mat'))
%         load(strcat('sim4_result_beam_gridonly_gs',num2str(gridsize),'.mat'))
%         load(strcat('sim4_result_random_gridonly_gs',num2str(gridsize),'.mat'))
%     end
    metrics_all(k,:,:) = metrics;
    metrics_beam_all(k,:,:) = metrics_beam;
    metrics_random_all(k,:,:) = metrics_random;

    if k<K-1
        load(strcat('sim4_result_seedsize5_gs',num2str(gridsize),'.mat'))
        load(strcat('sim4_result_random_seedsize5_gs',num2str(gridsize),'.mat'))
        load(strcat('sim4_result_beam_seedsize5_gs',num2str(gridsize),'.mat'))
    else

    end
    metrics_all_fix(k,:,:) = metrics;
    metrics_beam_all_fix(k,:,:) = metrics_beam;
    metrics_random_all_fix(k,:,:) = metrics_random;
end


% offsets = [-0.25,-0.05,0.15;-0.25,-0.05,0.15;-0.25,-0.05,0.15;-0.25,-0.05,0.15]';
% offsets_fix = [-0.15,0.05,0.25;-0.15,0.05,0.25;-0.25,-0.05,0.15;-0.25,-0.05,0.15]';

offsets = [-0.25,-0.05,0.15;-0.25,-0.05,0.15;-0.25,-0.05,0.15]';
offsets_fix = [-0.15,0.05,0.25;-0.15,0.05,0.25;-0.25,-0.05,0.15]';

colors = ['k','b', 'r'];
qt = [25,75];
tcl = tiledlayout(3,3);
tcl.TileSpacing = 'compact';
tcl.Padding = 'compact';
gridsname = {'17x17','9x9','5x5'};%{'17x17','9x9','5x5','5x5 /local'};%[0.05,0.1,0.2];
for i = 1:length(sample_factors)*length(factors)
nexttile(tcl)
    ax = multiple_merge_plot(gca, metrics_all(:,i,:), ...
        metrics_beam_all(:,i,:),metrics_random_all(:,i,:), ...
        gridsname,qt,colors,offsets,'-','o',{'greedy','beam','random'});
    ax = multiple_merge_plot(ax, metrics_all_fix(:,i,:), ...
        metrics_beam_all_fix(:,i,:),metrics_random_all_fix(:,i,:), ...
        gridsname,qt,colors,offsets_fix,'-.','^',{'greedy+seedsize5','beam+seedsize5','random+seedsize5'});
    min_white_margin(gca, 0.2, 0.1)

    if mod(i-1,length(sample_factors))+1==1
        ylabel(ax,strcat('\sigma = ',num2str(factors(ceil(i/length(factors))))))
    end
    if ceil(i/length(factors))==3
        xlabel(ax,strcat('\beta = ',num2str(sample_factors(mod(i-1,length(sample_factors))+1))));
    end
end
hL = legend(ax); 
% Move the legend to the right side of the figure
hL.Layout.Tile = 'East';
set(gcf, 'Position', [50 50 600 500]); %



figure
offsets_fix = [-0.2,0.0,0.2;-0.2,0.0,0.2;-0.2,0.0,0.2]';

colors = ['k','b', 'r'];
qt = [0,100];
tcl = tiledlayout(3,3);
tcl.TileSpacing = 'compact';
tcl.Padding = 'compact';
gridsname = {'17x17','9x9','5x5'};%{'17x17','9x9','5x5','5x5 /local'};%[0.05,0.1,0.2];
for i = 1:length(sample_factors)*length(factors)
nexttile(tcl)
    ax = multiple_merge_boxplot(gca, metrics_all(:,i,:), ...
        metrics_beam_all(:,i,:),metrics_random_all(:,i,:), ...
        gridsname,qt,colors,offsets,'-','o',{'greedy','beam','random'});
%     ax = multiple_merge_plot(ax, metrics_all_fix(:,i,:), ...
%         metrics_beam_all_fix(:,i,:),metrics_random_all_fix(:,i,:), ...
%         gridsname,qt,colors,offsets_fix,'-.','^',{'greedy+seedsize5','beam+seedsize5','random+seedsize5'});
    min_white_margin(gca, 0.2, 0.1)

    if mod(i-1,length(sample_factors))+1==1
        ylabel(ax,strcat('\sigma = ',num2str(factors(ceil(i/length(factors))))))
    end
    if ceil(i/length(factors))==3
        xlabel(ax,strcat('\beta = ',num2str(sample_factors(mod(i-1,length(sample_factors))+1))));
    end
end
% hL = legend(ax); 
% % Move the legend to the right side of the figure
% hL.Layout.Tile = 'East';
set(gcf, 'Position', [50 50 500 500]); %
saveas(gcf, 'merge_compare_boxplot.eps', 'epsc')


figure
offsets_fix = [-0.2,0.0,0.2;-0.2,0.0,0.2;-0.2,0.0,0.2]';

colors = ['k','b', 'r'];
qt = [0,100];
tcl = tiledlayout(3,3);
tcl.TileSpacing = 'compact';
tcl.Padding = 'compact';
gridsname = {'17x17','9x9','5x5'};%{'17x17','9x9','5x5','5x5 /local'};%[0.05,0.1,0.2];
for i = 1:length(sample_factors)*length(factors)
nexttile(tcl)
%     ax = multiple_merge_boxplot(gca, metrics_all(:,i,:), ...
%         metrics_beam_all(:,i,:),metrics_random_all(:,i,:), ...
%         gridsname,qt,colors,offsets,'-','o',{'greedy','beam','random'});
    ax = multiple_merge_boxplot(ax, metrics_all_fix(:,i,:), ...
        metrics_beam_all_fix(:,i,:),metrics_random_all_fix(:,i,:), ...
        gridsname,qt,colors,offsets_fix,'-.','^',{'greedy','beam','random'});
    min_white_margin(gca, 0.2, 0.1)

    if mod(i-1,length(sample_factors))+1==1
        ylabel(ax,strcat('\sigma = ',num2str(factors(ceil(i/length(factors))))))
    end
    if ceil(i/length(factors))==3
        xlabel(ax,strcat('\beta = ',num2str(sample_factors(mod(i-1,length(sample_factors))+1))));
    end
end
% hL = legend(ax); 
% % Move the legend to the right side of the figure
% hL.Layout.Tile = 'East';
set(gcf, 'Position', [50 50 500 500]); %
saveas(gcf, 'merge_compare_seed5_boxplot.eps', 'epsc')
