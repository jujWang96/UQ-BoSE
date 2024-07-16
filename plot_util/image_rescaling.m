function ax = image_rescaling(ax,inputname, add_xlabel, add_ylabel, keep_offset)
    
    function y = calc_sec(x,aboxwdt,boxmin,boxcen,darcsecperpix,offset,Axis)
        y = x-offset;
        if Axis==1
            y = -y/darcsecperpix;
        else
            y = y/darcsecperpix;
        end
        y = y+boxcen;
        y = y-boxmin;
        y = y/aboxwdt;
    end
box on
if strcmp('NGC2300',inputname) || strcmp('NGC2300GCR',inputname) ||strcmp('NGC2300centroid',inputname) 
    calc_ra = @(x) calc_sec(x,706.0,3062.0,3415.0,0.492,0,1);
    calc_dec = @(x) calc_sec(x,682.0,3489.0,3830.0,0.492,0,2);
elseif strcmp('NGC2300XMM',inputname) || strcmp('NGC2300XMMGCR',inputname)||strcmp('NGC2300XMMcentroid',inputname) 
    calc_ra = @(x) calc_sec(x,4800.0,23707.5,26107.5,0.05,5.38516,1);
    calc_dec = @(x) calc_sec(x,4800.0,25646.5,28046.5,0.05,-6.92139,2);
elseif strcmp('Arp299XMM',inputname) || strcmp('Arp299XMMGCR',inputname)||strcmp('Arp299XMMcentroid',inputname) 
    calc_ra = @(x) calc_sec(x,7040,21640.5,25160.5,0.05,2.80783,1);
    calc_dec = @(x) calc_sec(x,7044.0,24321.5,27843.5,0.05,-4.02374,2);

else 
   stop('filename is incorrect')
end
    axis([0 1 0 1 ])
    set(gcf, 'Position', [50 50 400 400]); %
    r = 1024;

    if strcmp('NGC2300',inputname) || strcmp('NGC2300XMM',inputname) || strcmp('Arp299XMM',inputname)
        set(gca,'XTick',[calc_ra(50),calc_ra(40),calc_ra(30),calc_ra(20),calc_ra(10),...
            calc_ra(0),calc_ra(-10),calc_ra(-20),calc_ra(-30),calc_ra(-40),calc_ra(-50)],...
            'YTick', [calc_dec(-50),calc_dec(-40),calc_dec(-30),calc_dec(-20),calc_dec(-10),...
            calc_dec(0),calc_dec(10),calc_dec(20),calc_dec(30),calc_dec(40),calc_dec(50)],...
            'Fontsize',18,'FontWeight','bold')
        xticklabels({'50','40','30','20','10','0','-10','-20','-30','-40','-50'})
        yticklabels({'-50','-40','-30','-20','-10','0','10','20','30','40','50'})
        xlim([calc_ra(60),calc_ra(-60)])
        ylim([calc_dec(-60),calc_dec(60)])
        
    elseif strcmp('NGC2300GCR',inputname) || strcmp('NGC2300XMMGCR',inputname) ||strcmp('Arp299GCR',inputname)
        set(gca,'XTick',[calc_ra(50),calc_ra(40),calc_ra(30),calc_ra(20),calc_ra(10),...
            calc_ra(0),calc_ra(-10),calc_ra(-20),calc_ra(-30),calc_ra(-40),calc_ra(-50)]*r,...
            'YTick', flip(r-[calc_dec(-50),calc_dec(-40),calc_dec(-30),calc_dec(-20),calc_dec(-10),...
            calc_dec(0),calc_dec(10),calc_dec(20),calc_dec(30),calc_dec(40),calc_dec(50)]*r),'Fontsize',18,'FontWeight','bold')
        set(gca,'TickDir','in')
        
        xticklabels({'50','40','30','20','10','0','-10','-20','-30','-40','-50'})
        yticklabels(flip({'-50','-40','-30','-20','-10','0','10','20','30','40','50'}))
        xlim([calc_ra(60),calc_ra(-60)]*r)
        ylim([1-calc_dec(60),1-calc_dec(-60)]*r)
        
    elseif strcmp('Arp299XMM',inputname) || strcmp('Arp299XMMGCR',inputname) ||strcmp('Arp299GCR',inputname)
        set(gca,'XTick',[calc_ra(50),calc_ra(40),calc_ra(30),calc_ra(20),calc_ra(10),...
            calc_ra(0),calc_ra(-10),calc_ra(-20),calc_ra(-30),calc_ra(-40),calc_ra(-50)]*r,...
            'YTick', flip(r-[calc_dec(-50),calc_dec(-40),calc_dec(-30),calc_dec(-20),calc_dec(-10),...
            calc_dec(0),calc_dec(10),calc_dec(20),calc_dec(30),calc_dec(40),calc_dec(50)]*r),'Fontsize',18,'FontWeight','bold')
        set(gca,'TickDir','in')
        
        xticklabels({'50','40','30','20','10','0','-10','-20','-30','-40','-50'})
        yticklabels(flip({'-50','-40','-30','-20','-10','0','10','20','30','40','50'}))
        xlim([calc_ra(60),calc_ra(-60)]*r)
        ylim([1-calc_dec(60),1-calc_dec(-60)]*r)
    elseif strcmp('NGC2300centroid',inputname) || strcmp('NGC2300XMMcentroid',inputname) ||strcmp('Arp299XMMcentroid',inputname)
       set(gca,'XTick',[calc_ra(50),calc_ra(40),calc_ra(30),calc_ra(20),calc_ra(10),...
            calc_ra(0),calc_ra(-10),calc_ra(-20),calc_ra(-30),calc_ra(-40),calc_ra(-50)],...
            'YTick', [calc_dec(-50),calc_dec(-40),calc_dec(-30),calc_dec(-20),calc_dec(-10),...
            calc_dec(0),calc_dec(10),calc_dec(20),calc_dec(30),calc_dec(40),calc_dec(50)],...
            'Fontsize',18,'FontWeight','bold')
        xticklabels({'50','40','30','20','10','0','-10','-20','-30','-40','-50'})
        yticklabels({'-50','-40','-30','-20','-10','0','10','20','30','40','50'})
        xlim([calc_ra(30),calc_ra(-30)])
        ylim([calc_dec(-30),calc_dec(30)])
        
    end
    if strcmp('NGC2300',inputname)
        if add_xlabel
            if keep_offset
                xlabel('\DeltaRA [arcsec] +7:32:16.0517')
            else
                xlabel('\DeltaRA [arcsec]')
            end
        else
            xlabel(' ')
        end
        if add_ylabel
            if keep_offset
                ylabel('\DeltaDec [arcsec] +85:42:39.101')
            else
                ylabel('\DeltaDec [arcsec]')
            end
        else
            ylabel(' ')
        end
        text(0.45, 0.65, 'NGC2300 (Chandra)','FontSize',18,'FontWeight','bold')
        set(gca,'XTickLabelRotation', 0,'linewidth',2)
        set(gcf,'renderer','Painters')    

    elseif strcmp('NGC2300XMM',inputname)
        if add_xlabel
            if keep_offset
                xlabel('\DeltaRA [arcsec] +7:32:20.8511')
            else
                xlabel('\DeltaRA [arcsec]')
            end
        else
            xlabel(' ')
        end
        if add_ylabel
            if keep_offset
                ylabel('\DeltaDec [arcsec] +85:42:32.186')
            else
                ylabel('\DeltaDec [arcsec]')
            end
        else
            ylabel(' ')
        end
        text(0.5, 0.75, 'NGC2300 (XMM)','FontSize',18,'FontWeight','bold')
        set(gca,'XTickLabelRotation', 0,'linewidth',2)
        set(gcf,'renderer','Painters')    

    elseif strcmp('Arp299XMM',inputname)
        if add_xlabel
            if keep_offset
                xlabel('\DeltaRA [arcsec] +11:28:32.3109')
            else
                xlabel('\DeltaRA [arcsec]')
            end
        else
            xlabel(' ')
        end
        if add_ylabel   
            if keep_offset
                ylabel('\DeltaDec [arcsec] +58:33:43.944')
            else
                ylabel('\DeltaDec [arcsec]')
            end
        else
            ylabel(' ')
        end
        text(0.52, 0.65, 'Arp299 (XMM)','FontSize',18,'FontWeight','bold')
        set(gca,'XTickLabelRotation', 0,'linewidth',2)
        set(gcf,'renderer','Painters')    

    elseif strcmp('NGC2300GCR',inputname)
        if add_xlabel
            if keep_offset
                xlabel('\DeltaRA [arcsec] +7:32:16.0517')
            else
                xlabel('\DeltaRA [arcsec]')
            end
        else
            xlabel(' ')
        end
        if add_ylabel
            if keep_offset
                ylabel('\DeltaDec [arcsec] +85:42:39.101')
            else
                ylabel('\DeltaDec [arcsec]')
            end
        else
            ylabel(' ')
        end
        text(0.45*r, (1-0.65)*r, 'NGC2300 (Chandra)','FontSize',18,'FontWeight','bold')
        set(gca,'XTickLabelRotation', 0,'linewidth',2)
        set(gcf,'renderer','Painters')    
        
    elseif strcmp('NGC2300XMMGCR',inputname)  
        if add_xlabel
            if keep_offset
                xlabel('\DeltaRA [arcsec] +7:32:20.8511')
            else
                xlabel('\DeltaRA [arcsec]')
            end
        else
            xlabel(' ')
        end
        if add_ylabel   
            if keep_offset
                ylabel('\DeltaDec [arcsec] +85:42:32.186')
            else
                ylabel('\DeltaDec [arcsec]')
            end
        else
            ylabel(' ')
        end
        text(0.51*r, (1-0.75)*r, 'NGC2300 (XMM)','FontSize',18,'FontWeight','bold')
        set(gca,'XTickLabelRotation', 0,'linewidth',2)
        set(gcf,'renderer','Painters') 
    elseif strcmp('Arp299XMMGCR',inputname)  
        if add_xlabel
            if keep_offset
                xlabel('\DeltaRA [arcsec] +11:28:32.3109')
            else
                xlabel('\DeltaRA [arcsec]')
            end
        else
            xlabel(' ')
        end
        if add_ylabel
            if keep_offset
                ylabel('\DeltaDec [arcsec] +58:33:43.944')
            else
                ylabel('\DeltaDec [arcsec]')
            end
        else
            ylabel(' ')
        end
        text(0.5*r, (1-0.65)*r, 'ARP299 (XMM)','FontSize',18,'FontWeight','bold')
        set(gca,'XTickLabelRotation', 0,'linewidth',2)
        set(gcf,'renderer','Painters')   
    elseif strcmp('NGC2300centroid',inputname)
        Ylm=ylim;                          % get x, y axis limits 
        Xlm=xlim;                          % so can position relative instead of absolute
        Xlb=mean(Xlm);                    % set horizontally at midpoint
        Ylb=0.99*Ylm(1);      
%         xlabel('\DeltaRA [arcsec] +7:32:16.0517', ...
%             'Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center')
        xlabel('\DeltaRA [arcsec]', ...
            'Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center')
        Xlb=0.98*Xlm(1);                    % set horizontally at midpoint
        Ylb = mean(Ylm); 
%         ylabel('\DeltaDec [arcsec] +85:42:39.101','Position',[Xlb Ylb], ...
%             'VerticalAlignment','baseline','HorizontalAlignment','center')
        ylabel('\DeltaDec [arcsec]','Position',[Xlb Ylb], ...
            'VerticalAlignment','baseline','HorizontalAlignment','center')
        text(0.47, 0.57, 'NGC2300 (Chandra)','FontSize',18,'FontWeight','bold')
        set(gca,'XTickLabelRotation', 0,'linewidth',2)
        set(gcf,'renderer','Painters')    
    elseif strcmp('NGC2300XMMcentroid',inputname)
        Ylm=ylim;                          % get x, y axis limits 
        Xlm=xlim;                          % so can position relative instead of absolute
        Xlb=mean(Xlm);                    % set horizontally at midpoint
        Ylb=0.99*Ylm(1);      
%         xlabel('\DeltaRA [arcsec] +7:32:20.8511', ...
%             'Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center')
        xlabel('\DeltaRA [arcsec]', ...
            'Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center')
        Xlb=0.98*Xlm(1);                    % set horizontally at midpoint
        Ylb = mean(Ylm); 
%         ylabel('\DeltaDec [arcsec] +85:42:32.186','Position',[Xlb Ylb], ...
%             'VerticalAlignment','baseline','HorizontalAlignment','center')
        ylabel('\DeltaDec [arcsec]','Position',[Xlb Ylb], ...
            'VerticalAlignment','baseline','HorizontalAlignment','center')
        text(0.51, 0.63, 'NGC2300 (XMM)','FontSize',18,'FontWeight','bold')
        set(gca,'XTickLabelRotation', 0,'linewidth',2) 
        set(gcf,'renderer','Painters')    

    elseif strcmp('Arp299XMMcentroid',inputname)
        Ylm=ylim;                          % get x, y axis limits 
        Xlm=xlim;                          % so can position relative instead of absolute
        Xlb=mean(Xlm);                    % set horizontally at midpoint
        Ylb=0.99*Ylm(1);      
%         xlabel('\DeltaRA [arcsec] +11:28:32.3109', ...
%             'Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center')
        xlabel('\DeltaRA [arcsec]', ...
            'Position',[Xlb Ylb],'VerticalAlignment','top','HorizontalAlignment','center')
        Xlb=0.98*Xlm(1);                    % set horizontally at midpoint
        Ylb = mean(Ylm); 
%         ylabel('\DeltaDec [arcsec] +58:33:43.944','Position',[Xlb Ylb], ...
%             'VerticalAlignment','baseline','HorizontalAlignment','center')
        ylabel('\DeltaDec [arcsec]','Position',[Xlb Ylb], ...
            'VerticalAlignment','baseline','HorizontalAlignment','center')
        text(0.5, 0.58, 'ARP200 (XMM)','FontSize',18,'FontWeight','bold')
        set(gca,'XTickLabelRotation', 0,'linewidth',2)    
        set(gcf,'renderer','Painters')    

    end
%     % Define the desired margin sizes in inches
%     leftMargin = 20;   % Left margin
%     rightMargin = 1.0;  % Right margin
%     topMargin = 0.5;    % Top margin
%     bottomMargin = 20; % Bottom margin
%     
%     % Get the current figure position and size
%     figPos = get(gcf, 'Position');
%     
%     % Calculate the new figure position and size to set the margins
%     newWidth = figPos(3) + leftMargin + rightMargin;
%     newHeight = figPos(4) + topMargin + bottomMargin;
%     newLeft = figPos(1) - leftMargin;
%     newBottom = figPos(2) - bottomMargin;
%     
%     % Set the new figure position and size
%     set(gcf, 'Position', [newLeft, newBottom, newWidth, newHeight]);
%     
%         
    %axis equal

end