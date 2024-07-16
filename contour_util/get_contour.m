function [contour,N] = get_contour(n, DT,selected,adj_mat,invalid,cell_area)

% Get the contour of segmentations 
%Input: 
%
% n: number of points 
% DT: a Delaunay triangulation object
% index: the index of the segmentation with samllest BIC/or forced to
% choose
% sets_all: all segmentations
% adj_mat: adjacency matrix 
% invalid: cells with NaN log intensity and isolated cells
%
%
%
%
%
%Output:
% contour: contour object 
% N: number of segmented region 

[V, R] = voronoiDiagram(DT);
vx = V(:,1);
vy = V(:,2);
regionId = find(~cellfun(@isempty,selected));
if length(regionId)==1
    contour = {};
    N = 1;
    return
end
allV = 1:n;
contourV = {};
connectList = [];
%cellnum-by-2 matrix, 1 covers 2
coverList = [];
background=[];
point_on_bound = [];
for c = 1:n
    if isnan(cell_area(c))
        point_on_bound = [point_on_bound c];
    end
end

%check if a segment touches the domain. If so consider it as background. 
for seg = regionId
    if any(reshape(adj_mat(point_on_bound,selected{seg}),1,[]))
        background = [background seg];
    end
end
%disp('background check')
for i = invalid
    adj_mat(i,:) = false;
    adj_mat(:,i) = false;
    
end
N = 0;

for seg = regionId
        contourcurr = [];
        N = N+1;
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
        contourV{seg} = [vx(contourcurr(:,1)),vx(contourcurr(:,2)),vy(contourcurr(:,1)),vy(contourcurr(:,2))];
    
end
%check the connectivity and coverage
for seg1 = regionId
    for seg2 = regionId
        neighbourLine = length(intersect(contourV{seg1},contourV{seg2},'rows'));
        if seg2>seg1 && neighbourLine >1
            connectList = [connectList;[seg1, seg2]];
           
        end
    end

end
%calculate area and intensity 
segArea = [];
segIntensity = [];
for seg = regionId
    areaList = cell_area(selected{seg});
    areaListClean = areaList(~isnan(areaList));
    segArea = [segArea, sum(areaListClean)];
    segIntensity = [segIntensity, sum(~isnan(areaList))/sum(areaListClean)];
end

contour.contourV = contourV;
contour.regionId = regionId;
contour.background = background;
contour.connectList = connectList;
contour.segArea = segArea;
contour.segIntensity = segIntensity;
end
