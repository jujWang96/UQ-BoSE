function [area_seg,flux_seg,ratio_seg,x_seq,y_seq] = get_area_flux(X,selected,cell_area,sortflux)
area_seg = [];
flux_seg = [];
x_seq = [];
y_seq = [];
for ii = 1:length(selected)
    if ~isempty(selected{ii})
        area_seg = [area_seg, sum(cell_area(selected{ii}))];
        flux_seg = [flux_seg, length(selected{ii})];
        x_seq = [x_seq, mean(X((selected{ii}),1)) ];
        y_seq = [y_seq, mean(X((selected{ii}),2)) ];

    end  
end
ratio_seg = flux_seg./area_seg;

if sortflux
    [ratio_seg,I] = sort(ratio_seg);
    area_seg = area_seg(I);
    flux_seg = flux_seg(I);  
    x_seq = x_seq(I);
    y_seq = y_seq(I);

end

