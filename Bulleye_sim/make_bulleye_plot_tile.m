close all
clear
addpath(genpath('~/UQ-BoSE/Bulleye_sim/output'))
SNRs = [3,10,30];
linetypes = {'--','-.','-'};
linecolors = {'k','r','b'};
tcl = tiledlayout(3,3);
tcl.TileSpacing = 'compact';
tcl.Padding = 'compact';
displacenames = {'\mu_B = 100','\mu_B = 500','\mu_B = 2500'};
same_contrast = [1,4,7];
contrasts = [1.4,2,2.8];
lambdas = [100, 500,2500];
csvFileName = '~/UQ-BoSE/Bulleye_sim/output/normality.csv';
headers = {'SNR','contrast', 'lambda', 'valid_num','p-val','<=0.05'};
fid = fopen(csvFileName, 'w');
fprintf(fid, '%s,', headers{1:end-1});  % Write headers except the last one
fprintf(fid, '%s\n', headers{end});     % Write the last header with a newline
fclose(fid);

for i = 1:length(SNRs)
    addtitle = false;
    addxlabel = false;

    SNR = SNRs(i);
    load(strcat('sim_bulleye_result_SNR',num2str(SNR),'.mat'))
    rads = cell(1,9);
    for j = 1:9
        valid = find(rad(j,:)>0);
        lambda =lambdas(ceil(j/3));
        contrast = contrasts(mod(j-1,3)+1);
        %figure
        %ax = single_bulleye_plot(ax,rad(i,valid)','-');
        %figure
        rads{j} = rad(j,valid)';
        %boxplot(metrics(i,valid))
        disp(strcat("number of valid results:","SNR:",num2str(SNR),"; lambda:",num2str(lambda),"; constrast:",num2str(contrast)))
        disp(length(valid))
        disp(strcat("normality test result:","SNR:",num2str(SNR),"; lambda:",num2str(lambda),"; constrast:",num2str(contrast)))
        normalitytest(rad(j,valid))
        rad_std = normalize(rad(j,valid));
        disp(strcat("normality test result:","SNR:",num2str(SNR),"; lambda:",num2str(lambda),"; constrast:",num2str(contrast)))
        [ks,pval] = kstest(rad_std);
        pval
        f = figure;
        sgtitle(strcat("SNR:",num2str(SNR),"; lambda:",num2str(lambda),"; constrast:",num2str(contrast),"; normality test pval",num2str(pval)),'fontsize',15)

        f.Position = [100 100 600 250];
        subplot(1,2,1)
        %[f,xi] = ksdensity(rad(j,valid),'Support','positive','BoundaryCorrection','reflection','Bandwidth',0.01);
        %plot(xi,f);
        histfit(rad(j,valid),20,'kernel')
        subplot(1,2,2)
        qqplot(rad_std)
        saveas(gcf,strcat("~/UQ-BoSE/Bulleye_sim/output/normality-SNR",num2str(SNR),"lambda",num2str(lambda),"constrast",num2str(contrast*10)),'epsc')
        dataRow = [SNR, contrast,lambda,length(valid), pval, ks];
        dlmwrite(csvFileName, dataRow, '-append', 'precision', 6, 'delimiter', ',');

    end
if i==1
    addtitle = true;
end
if i==length(SNRs)
    addxlabel = true;
end
nexttile(tcl)
multiple_bulleye_plot(gca,rads(same_contrast),linetypes,linecolors,displacenames,addtitle,strcat('contrast = ',num2str(contrasts(1))),addxlabel,'|G_1|');
ylabel(strcat('SNR = ',num2str(SNR)))
nexttile(tcl)
multiple_bulleye_plot(gca,rads(same_contrast+1),linetypes,linecolors,displacenames,addtitle,strcat('contrast = ',num2str(contrasts(2))),addxlabel,'|G_1|');
nexttile(tcl)
ax = multiple_bulleye_plot(gca,rads(same_contrast+2),linetypes,linecolors,displacenames,addtitle,strcat('contrast = ',num2str(contrasts(3))),addxlabel,'|G_1|');
end
% hL = legend(ax); 
% 
% hL.Layout.Tile = 'East';

