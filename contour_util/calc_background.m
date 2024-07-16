function num = calc_background(selected,adj_mat,cell_area)

point_on_bound = [];
regionId = find(~cellfun(@isempty,selected));

for c = 1:length(cell_area)
    if isnan(cell_area(c))
        point_on_bound = [point_on_bound c];
    end
end
num=0;
for seg = regionId
    if any(reshape(adj_mat(point_on_bound,selected{seg}),1,[]))
        num = num+1;
    end
end
end