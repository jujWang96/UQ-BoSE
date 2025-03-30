clear
close all

parpath = '~/src/gsrg/';
plotpath = 'coverage_sim/plots/';
cd(strcat(parpath, 'coverage_sim'))

resultpath = 'coverage_sim/results/';
addpath(genpath(strcat(parpath,'Astro_sim')))
addpath(genpath(strcat(parpath,'G-SRG')))
addpath(genpath(strcat(parpath,'plot_util')))
addpath(genpath(strcat(parpath,'contour_util')))
addpath(genpath(strcat(parpath,'FD_util')))
addpath(genpath(strcat(parpath,'boot_util')))
addpath(genpath(strcat(parpath,'gcr_util')))


gridsize = 0.1;
seedsize = 5;
rep_itr = 1000;

grid_only = true;
force_merge = true;
pair_search = true;

merge_method = "random";
factor = 30;
sample_factor = 1;

shapenames = ["spiral"];%["ellipse","bridge","spiral"];
methods = ["BIC","true_FD_num"]; 
label_p = {texlabel('P_{LC}'),texlabel('P_G'),texlabel('P_l')};
label_psim = {texlabel('P^{BIC}_e'),texlabel('P_e')};
linedot_p = ["-","-.","--"];
linedot_psim = ["o-","*-"];
for shapename = shapenames
    p_sims_all = {};
    i = 0;
    for method = methods
        i = i+1;
        selection_method = method;
        boot_method = method;
        clear p_sims
        load(strcat(parpath,resultpath,"coverage_prob_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor)  , "_gridsize", num2str(gridsize), ...
                    "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search),...
                    "_selection_method",selection_method, "_boot_method" ,boot_method,'.mat'))
        p_sims_all{i} = p_sims;
    end
    ps_all{1} = p_LCs;
    ps_all{2} = p_Gs;
    ps_all{3} = p;
    range = [0.5,3];
    
    make_plot_coverage_prob(range,Beta,true,ps_all,p_sims_all,label_p, label_psim, linedot_p, linedot_psim,'k','k')
    saveas(gca,strcat(parpath,plotpath,"coverage_prob_log_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor)  , "_gridsize", num2str(gridsize), ...
                    "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search)),"epsc")
    
    make_plot_coverage_prob(range,Beta,false,ps_all,p_sims_all,label_p, label_psim, linedot_p, linedot_psim,'k','k')
    saveas(gca,strcat(parpath,plotpath,"coverage_prob_",shapename,"_sample_factor", num2str(sample_factor),"_factor", num2str(factor)  , "_gridsize", num2str(gridsize), ...
                    "_gridonly",string(grid_only),"_merge",merge_method,"_forcemerge",string(force_merge),"_pairsearch",string(pair_search)),"epsc")
end
GCRs = GCRs(~invalid_index);
origin_coef_set = origin_coef_set(:,~invalid_index);
zero_rows = all(covers == 0, 2);
zero_row_indices = find(zero_rows);
idx = 259;
figure
GCR= GCRs{idx};


plot(GCR{4})
hold on 
[invx,invy] = iFD(origin_coef_set(:,idx),true_FD_num);
plot(invx,invy,'r')

[invx,invy] = iFD(origin_coef_set(:,idx));
plot(invx,invy,'k')

plot(V(:,1),V(:,2),'k')
axis equal

% plot the distribution of fourier coefficient 
[rep,D] = size(para);
for d = 1:D
    if d <= ceil(D/2)
        fd = d-ceil(true_FD_num/2);
        fd_text = strcat("Real",num2str(fd));
    else
        fd = d-ceil(D/2)-ceil(true_FD_num/2);
        if fd >= 1
            fd = fd+1;
        end
        fd_text = strcat("Imag",num2str(fd));

    end
    check_normality(real(para(:,d)),d,fd_text)
    saveas(gcf,strcat(parpath,plotpath,shapenames,"normality-d",num2str(d),"-fd",fd_text),'epsc')

end

filename = strcat(parpath,resultpath,'coef.csv');

% Write the matrix to the CSV file
writematrix(para, filename);
