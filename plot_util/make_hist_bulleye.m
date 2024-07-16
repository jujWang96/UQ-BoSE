function [] = make_hist_bulleye(SNR,c,base,valid_only,grid,path)
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.08 0.02], [0.05 0.02]);
idx=0;
if nargin<6
    path = '~/Desktop/QE/simple_simulation/Bulleye/';
end

for i = 1:length(c)
    for j = 1:length(base)     
        idx = idx+1;
        if nargin==5
            load(strcat(path,'bulleye_area_base', num2str(base(j)),'c',num2str(c(i)),'SNR',num2str(SNR),'grid',num2str(grid(j)),'m4.mat'))
        else
            load(strcat(path,'bulleye_area_base', num2str(base(j)),'c',num2str(c(i)),'SNR',num2str(SNR),'m4.mat'))

        end
        if valid_only 
            rad = rad(rad>0);
            valid_prop = length(rad)/500;
        end
        subplot(length(c),length(base),idx)
        subplot_handle = subplot(length(c),length(base),idx);
        subplot_pos = get(subplot_handle, 'position');
        histogram(rad,'BinEdges', linspace(-0.01,0.5,51),'Normalization','probability');
        xlim([0,0.3])
        %ylim([0,1])
        if i==length(c)
            h=xlabel('|G_1|');
            set(h, 'FontSize', 15) ;

        end
    add_annotation([0.7 0.75 0.08 0.08], subplot_pos, ['c=', num2str(c(i))])
    add_annotation([0.7 0.6 0.08 0.08], subplot_pos, ['\alpha=', num2str(base(j)/100)])
    if valid_only
        add_annotation([0 0.75 0.08 0.08], subplot_pos, [ num2str(round(valid_prop*100,2)),'%'])
    end

    set(gca,'FontSize',12)

    end



end