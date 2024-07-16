function [contour,selected] = merge_segments(contour,selected,seg_to_merge)
%update segmentation and contour after merging all candidate segments
%input:
% seg_to_merge:candidate segments to be merged
target_seg = seg_to_merge(1);
origin_regionId = contour.regionId ;

%merge the segments and update contourV
for i = 2:numel(seg_to_merge)
   
    contour.contourV{target_seg} = [setdiff(contour.contourV{target_seg},...
    contour.contourV{seg_to_merge(i)},'rows');setdiff(contour.contourV{seg_to_merge(i)},contour.contourV{target_seg},'rows')];
    selected{target_seg}=[selected{target_seg} selected{seg_to_merge(i)}];
    selected{seg_to_merge(i)}=[];
end
%update the contour.regionId
contour.regionId = find(~cellfun(@isempty,selected));
%update the contour.background
if ~isempty(intersect(seg_to_merge,contour.background))
    contour.background = setdiff(contour.background,seg_to_merge);
    contour.background = [contour.background,target_seg];
end
end
