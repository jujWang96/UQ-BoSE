function res = check_same_seg(seg1,seg2,num)
%input:
% seg1,seg2: 1 by num cell array representing two different segmentations,
% each cell represents a connected subsets
%output:
% True if the two segmentations are identical 


res = true;
for i = 1:num
    if isempty(seg1{i}) && isempty(seg2{i})
        continue
    elseif isempty(setxor(seg1{i},seg2{i}))
        continue
    else
        res = false;
        break
    end
end