function [] = make_plot_coverage_prob(range,Beta,logy,ps_all,p_sims_all,label_p, label_psim, linedot_p, linedot_psim,colorp,colorpe)
figure
if logy
    for ii =1:length(ps_all)
        semilogy(Beta,1-ps_all{ii},linedot_p(ii),'Color',colorp);
        hold on
    end
    for ii = 1:length(p_sims_all)
        semilogy(Beta,1-p_sims_all{ii},linedot_psim(ii),'Color',colorpe);
        hold on

    end
else
    for ii =1:length(ps_all)
        plot(Beta,1-ps_all{ii},linedot_p(ii),'Color',colorp);
        hold on

    end
    for ii = 1:length(p_sims_all)
        plot(Beta,1-p_sims_all{ii},linedot_psim(ii),'Color',colorpe);
        hold on

    end

end

xlim(range)
labels = [label_p,label_psim];
legend(labels,'Location','southwest','FontSize',15)
xlabel(texlabel('gamma'),'FontSize', 20)
ylabel('Error Probability','FontSize', 15)

set(gcf,'Position',[100 75 500 375])
set(gca,'FontSize',15)





% function [] = make_plot_coverage_prob(range,Beta,p,p_LCs,p_Gs,p_sim_BICs,p_sim_trues,logy,colorp,colorpe)
% 
% if logy
%     plot(Beta,1-p_LCs,'-','Color',colorp);
%     hold on     
% else
%     %log scale y for more detailed tail behavior
%     semilogy(Beta,1-p_LCs,'-','Color',colorp);
%     hold on 
% end
% plot(Beta,1-p_Gs,'-.','Color',colorp);
% 
% plot(Beta,1-p,'--','Color','k')
% plot(Beta,1-p_sim_BICs,'o-','Color',colorpe);
% plot(Beta,1-p_sim_trues,'*-','Color',colorpe)
% xlim(range)
% legend({texlabel('p_{LC}'),texlabel('p_G'),texlabel('p_l'),texlabel('p^{BIC}_e'),texlabel('p_e')},'Location','southwest','FontSize',12)
% xlabel(texlabel('gamma'),'FontSize', 20)
% ylabel('Error Probability','FontSize', 15)
% 
% set(gcf,'Position',[100 75 500 375])
% set(gca,'FontSize',15)



