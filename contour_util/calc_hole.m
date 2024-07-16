function num_of_hole  = calc_hole(selected,V,R,adj_mat,invalid,cell_area)
vx = V(:,1);
vy = V(:,2);
num_of_hole=0;
n = length(vx);
regionId = find(~cellfun(@isempty,selected));
if length(regionId)==1
    return
end
background = [];
point_on_bound = [];
regionId = find(~cellfun(@isempty,selected));

for c = 1:length(cell_area)
    if isnan(cell_area(c))
        point_on_bound = [point_on_bound c];
    end
end
for seg = regionId
    if any(reshape(adj_mat(point_on_bound,selected{seg}),1,[]))
        background = [background,seg];
    end
end
allV = 1:n;
for i = invalid
    adj_mat(i,:) = false;
    adj_mat(:,i) = false;
    
end

for seg = regionId
    clear cx cy
        contourcurr = [];
        R_i = selected{seg};
 
        for j = R_i
            for k = R_i
                adj_mat(j,k) = false;
                adj_mat(k,j) = false;                      
            end
            adjV = allV(adj_mat(j,:));
  
            for l = adjV
                contourcurr = [contourcurr; intersect(R{j}, R{l})];
                adj_mat(j,l) = false;
                %adj_mat(l,j) = false;
            end
        end
        x1 = vx(contourcurr(:,1));
        x2 = vx(contourcurr(:,2));
        y1 = vy(contourcurr(:,1));
        y2 = vy(contourcurr(:,2));
    cx(1) = x1(1);
    cx(2) = x2(1);
    
    cy(1) = y1(1);
    cy(2) = y2(1);
    
    x1(1) = nan;
    x2(1) = nan;
    y1(1) = nan;
    y2(1) = nan;
    i=2;
    n_seg = length(x1);
    while i<n_seg
    %for i = 2:(n-1)
    try
        j = find(x1==cx(i) & y1==cy(i));
        if ~isempty(j)
            cx(i+1) = x2(j);
            cy(i+1) = y2(j);
        else
            j = find(x2==cx(i) & y2==cy(i));
            cx(i+1) = x1(j);
            cy(i+1) = y1(j);
        end
        i = i+1;
         x1(j) = nan;
         x2(j) = nan;
         y1(j) = nan;
         y2(j) = nan;
    catch
        if ~ any(background==seg)
            num_of_hole = num_of_hole+1;
        end
        j = find(~isnan(x1), 1);
        cx(i+1) = x1(j);
        cy(i+1) = y1(j);
        cx(i+2) = x2(j);
        cy(i+2) = y2(j);
        i = i+2;
         x1(j) = nan;
         x2(j) = nan;
         y1(j) = nan;
         y2(j) = nan;
    end
    end
    
end