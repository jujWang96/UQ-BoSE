
close all;
addpath(genpath('/Users/joycewang/Desktop/gsrg/Astro_sim/sim'))
addpath(genpath('/Users/joycewang/Desktop/gsrg/plot_util'))

T  = 500;
g1 = [ones(1500, 1); 2*ones(1500, 1); 3*ones(1500, 1)];
g2 = [ones(T, 3); 2*ones(T, 3); 3*ones(T, 3)];
g2 = g2(:);
reduction = 0.02;
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.04], [0.1 0.02], [0.08 0.02]);

colors = ["black","blue","red"];



% plot for sim4 for median
offset = [-0.2,-0.05,0.15];
factors = [0.5,1,2];
gridsizes = [0.05,0.1,0.2];
K = length(gridsizes);
qt = [25,75];
figure;
ylims = [0,0.6;0.2,0.6;0.4,0.8]; 
for i = 1:3
    if i == length(factors)
        add_xlabel = true;
    else
        add_xlabel = false;
    end
    for j = 0:2
        subplot(3,3,i+j*3)

        hold on
        medval = zeros(3,K);
        for  k = 1:K
            gridsize = gridsizes(k);
            if gridsize==0
                load(strcat('sim4_all_seeds_result.mat'))
                load(strcat('sim4_all_seeds_result_beam.mat'))
                load(strcat('sim4_all_seeds_result_random.mat'))


            else
                load(strcat('sim4_result_gs',num2str(gridsize),'.mat'))
                load(strcat('sim4_result_random_gs',num2str(gridsize),'.mat'))
                load(strcat('sim4_result_beam_gs',num2str(gridsize),'.mat'))
            end
            medval(1,k) = median(metrics(i+j*3,:)); 
            medval(2,k) = median(metrics_beam(i+j*3,:)); 
            medval(3,k) = median(metrics_random(i+j*3,:)); 
         
            errorbar(k+offset(1),medval(1,k),prctile(metrics(i+j*3,:),qt(1))-medval(1,k),prctile(metrics(i+j*3,:),qt(2))-medval(1,k),colors(1))
            scatter(k+offset(1),medval(1,k),colors(1))
            errorbar(k+offset(2),medval(2,k),prctile(metrics_beam(i+j*3,:),qt(1))-medval(2,k),prctile(metrics_beam(i+j*3,:),qt(2))-medval(2,k),colors(2))
            scatter(k+offset(2),medval(2,k),colors(2))
            errorbar(k+offset(3),medval(3,k),prctile(metrics_random(i+j*3,:),qt(1))-medval(3,k),prctile(metrics_random(i+j*3,:),qt(2))-medval(3,k),colors(3))
            scatter(k+offset(3),medval(3,k),colors(3))
        end
        plot((1:K)+offset(1),medval(1,:),colors(1))
        plot((1:K)+offset(2),medval(2,:),colors(2))
        plot((1:K)+offset(3),medval(3,:),colors(3))
        xticks(1:K)
        %ylim(ylims(i,:))
        xticklabels(gridsizes)
        xlim([1+offset(1)*2,K+offset(k)*2])
        if j==2 
            %xlable(strcat('\beta = ',num2str(sample_factors(i))))
            xlabel(strcat('\beta =  ',num2str(sample_factors(i))))
        end

        if i==1 
            %xlable(strcat('\beta = ',num2str(sample_factors(i))))
            ylabel(strcat('\sigma =  ',num2str(factors(j+1))))
        end
        box on

    end

end

offset = [-0.15,0.0,0.2;-0.15,0.0,0.2;-0.2,-0.05,0.15]';

for i = 1:3
    if i == length(factors)
        add_xlabel = true;
    else
        add_xlabel = false;
    end
    for j = 0:2

        subplot(3,3,i+j*3)
        %min_white_margin(gca, -reduction/2, reduction)

        hold on
        medval = zeros(3,K);
        for  k = 1:K
            gridsize = gridsizes(k);
            if gridsize==0
                load(strcat('sim4_all_seeds_result.mat'))
                %load(strcat('sim4_all_seeds_result_beam.mat'))
                load(strcat('sim4_all_seeds_result_random.mat'))


            else
                if gridsize==0.2
                    load(strcat('sim4_result_gs',num2str(gridsize),'.mat'))
                    load(strcat('sim4_result_random_gs',num2str(gridsize),'.mat'))
                    load(strcat('sim4_result_beam_gs',num2str(gridsize),'.mat'))

                else
                    load(strcat('sim4_result_seedsize5_gs',num2str(gridsize),'.mat'))
                    load(strcat('sim4_result_random_seedsize5_gs',num2str(gridsize),'.mat'))
                    load(strcat('sim4_result_beam_seedsize5_gs',num2str(gridsize),'.mat'))
                end
            end
            medval(1,k) = median(metrics(i+j*3,:)); 
            medval(2,k) = median(metrics_beam(i+j*3,:)); 
            medval(3,k) = median(metrics_random(i+j*3,:)); 
         
            errorbar(k+offset(1,k),medval(1,k),prctile(metrics(i+j*3,:),qt(1))-medval(1,k),prctile(metrics(i+j*3,:),qt(2))-medval(1,k),colors(1))
            scatter(k+offset(1,k),medval(1,k),colors(1),"^")
            errorbar(k+offset(2,k),medval(2,k),prctile(metrics_beam(i+j*3,:),qt(1))-medval(2,k),prctile(metrics_beam(i+j*3,:),qt(2))-medval(2,k),colors(2))
            scatter(k+offset(2,k),medval(2,k),colors(2),"^")
            errorbar(k+offset(3,k),medval(3,k),prctile(metrics_random(i+j*3,:),qt(1))-medval(3,k),prctile(metrics_random(i+j*3,:),qt(2))-medval(3,k),colors(3))
            scatter(k+offset(3,k),medval(3,k),colors(3),"^")
        end
        plot((1:K)+offset(1,:),medval(1,:),'LineStyle','-.','Color',colors(1))
        %plot((1:K)+offset(2),medval(2,:)',colors(2))
        plot((1:K)+offset(3,:),medval(3,:),'LineStyle','-.','Color',colors(3))
        xticks(1:K)
        %ylim(ylims(i,:))
        xticklabels(gridsizes)
        xlim([1+offset(1,1)*2,K+offset(k,1)*2])
        box on

        

    end
end

set(gcf, 'Position', [50 50 600 500]); %
