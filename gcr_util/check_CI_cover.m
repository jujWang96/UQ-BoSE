function output = check_CI_cover(V,BW)
%check whether the true shape lies entirely inside the confidence region
%Input:
%BM:binary mask of the global confidence region 
%V: L-by-2 matrix, X-Y coordinates of true shape contour
%
%BW = rot90(BW);
%BW = flip(rot90(flip(rot90(BW)),3),2);
%BW = flip(BW);
dim = size(BW);
r = dim(1);
%BW_flip = rot90(flip(BW,2)); %flip x and y coordinates


if any(BW(round(V(:,2)*r),round(V(:,1)*r))==0)
    output = 0;
else 
    output = 1;
end



