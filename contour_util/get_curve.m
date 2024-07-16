function [cx,cy] = get_curve(contour,clockwise)
%
% Get the closed curve of the segmentation in traversy order 
%Input:
% contour: line segments of the contour 
% clockwise: logic argument, if true then return the curve in clockwise
% manner
%
%Output: 
% cx
% cy
    if isempty(contour)
        cx = [];
        cy = [];
        curve = [];
        return
    end
    x1 = contour(:,1);
    
    x2 = contour(:,2);
    y1 = contour(:,3);
    y2 = contour(:,4);
    n = length(x1);
    cx = zeros(n,1);
    cy = zeros(n,1);

    cx(1) = x1(1);
    cx(2) = x2(1);
    
    cy(1) = y1(1);
    cy(2) = y2(1);
    
    x1(1) = nan;
    x2(1) = nan;
    y1(1) = nan;
    y2(1) = nan;
    for i = 2:(n-1)
%         [j,k] = find([x1, x2]== cx(i));
%         if length(j)~=1
%             [j2,k2] = find([y1, y2]== cy(i));
%             j = intersect(j, j2);
%             k = intersect(k, k2);
%         end
%         
%         if k == 1
%             cx(i+1) = x2(j);
%             cy(i+1) = y2(j);
%         else
%             cx(i+1) = x1(j);
%             cy(i+1) = y1(j);
%         end
%         x1(j) = nan;
%         x2(j) = nan;
%         y1(j) = nan;
%         y2(j) = nan;
        j = find(x1==cx(i) & y1==cy(i));
        if ~isempty(j)
            cx(i+1) = x2(j);
            cy(i+1) = y2(j);
        else
            j = find(x2==cx(i) & y2==cy(i));
            cx(i+1) = x1(j);
            cy(i+1) = y1(j);
        end
         x1(j) = nan;
         x2(j) = nan;
         y1(j) = nan;
         y2(j) = nan;
        
    end
    % return the curve in a clockwise manner 
    if clockwise==true
        if ispolycw(cx, cy)
            curve = [cx,cy];

            return;
        else 
            cx = flip(cx);
            cy = flip(cy);
               curve = [cx,cy];

            return;
        end
    % return the curve in counter clockwise manner    
    else 
        if ispolycw(cx, cy)
            cx = flip(cx);
            cy = flip(cy);
               curve = [cx,cy];

            return;
        else 
               curve = [cx,cy];

            return;
        end
        
    end
    
end




