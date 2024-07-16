function res = check_same_seg_restrict(seg1,seg2)
%input:
% seg1,seg2: 1 by num cell array representing two different segmentations,
% each cell represents a connected subsets
%output:
% True if the two segmentations are identical 


% res = true;
% nonempty1 = find(~cellfun(@isempty,seg1));
% nonempty2 = find(~cellfun(@isempty,seg2));
% %check nonempty cell first
% if ~all(nonempty1==nonempty2)
%     res = false;
% %check each single nonempty cell    
% else
%     for i = 1:length(nonempty1)
%         if all(seg1{nonempty1(i)}==seg2{nonempty2(i)})
%             continue
%         else
%             res = false;
%             break
%         end
%     end
% 
% end
nonempty1 = seg1(~cellfun(@isempty,seg1));
nonempty2 = seg2(~cellfun(@isempty,seg2));
%check nonempty cell first
if ~isequal(nonempty1,nonempty2)
    res = false;
else
    res = true;
end
