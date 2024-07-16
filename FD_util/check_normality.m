function check_normality(x,coef_parameter, fd_text)
    x_std = normalize(x);
    %disp(strcat("normality test result:","SNR:",num2str(SNR),"; lambda:",num2str(lambda),"; constrast:",num2str(contrast)))
    [ks,pval] = kstest(x_std);
    isnormal = (pval>0.05);
    %pval
    f = figure;
    sgtitle(strcat("coef: ",num2str(coef_parameter),"; FD: ",fd_text,"; KS test pval: ",num2str(pval),"; is normal:",num2str(isnormal)),'fontsize',15)

    f.Position = [100 100 600 250];
    subplot(1,2,1)
    %[f,xi] = ksdensity(rad(j,valid),'Support','positive','BoundaryCorrection','reflection','Bandwidth',0.01);
    %plot(xi,f);
    histfit(x,20,'kernel')
    subplot(1,2,2)
    qqplot(x_std)
