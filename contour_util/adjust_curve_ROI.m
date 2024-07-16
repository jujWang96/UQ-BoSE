function [cxs,cys] = adjust_curve_ROI(cxs,cys, ROI,clockwise)
%adjust the curve with a predefined ROI (polyshape)
%the curve should be restrained inside the ROI
curve_num = length(cxs);

for i = 1:curve_num
    poly_shape = polyshape(cxs{i}, cys{i});
    intersection_poly = intersect(poly_shape, ROI);
    [bx,by] = boundary(intersection_poly);
    if isempty(bx)
        cxs{i} = [];
        cys{i} = [];
        continue;
    end
    if clockwise
        [cx,cy] = poly2cw(bx,by);
    else
        [cx,cy] = poly2ccw(bx,by);
    end
    cxs{i} = cx;
    cys{i} = cy;

end
