function [target_x, target_y] = extract_curve(region_in_order, seg,n, DT,adj_mat,invalid,cell_area)
% Extract curve of the segmentation boundary 
% regionInOrder: specify the order of region in sequence from inner most to
% outer most so that do iteratively merge to extract boundaries
% seg: segmentation
% seg_countour: contour of segmentation 
% target_x, target_y: the x and y cordinates of each closed curve.
target_x = {};
target_y = {};

for ii = 1:(length(region_in_order)-1)
    temp_seg = seg;
    reg_in = region_in_order(ii);
    reg_out = region_in_order(ii+1);
    for jj = 1:length(region_in_order)
        if jj < ii 
            temp_seg{reg_in} = [temp_seg{reg_in}, temp_seg{region_in_order(jj)}];
            temp_seg{region_in_order(jj)} = [];
        end
        if jj > ii+1
            temp_seg{reg_out} = [temp_seg{reg_out}, temp_seg{region_in_order(jj)}];
            temp_seg{region_in_order(jj)} = [];
        end
    end
    temp_contour = get_contour(n, DT,temp_seg,adj_mat,invalid,cell_area);
    temp_contourV = temp_contour.contourV;
    temp_region_in = setdiff(temp_contour.regionId, temp_contour.background);
    [target_x{ii},target_y{ii}] = get_curve(temp_contourV{temp_region_in},false);
end

