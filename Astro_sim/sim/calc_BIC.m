function val = calc_BIC(X,contour,selected,m)
%calculate the BIC value of the selected segmentation 
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
adj_mat = get_adj_mat( E, n );
[invalid, valid] = get_invalid_cells(cell_log_intensity, adj_mat, n);

log_like = -sum(log(1:length(X)))-length(X)+length(invalid);
for i = reshape(contour.regionId,1,[])
    log_like = log_like + numel(selected{i})*log(numel(selected{i})/sum(cell_area(selected{i})));
end
val = -2*log_like+m*(length(contour.regionId)-1)*log(length(X));



