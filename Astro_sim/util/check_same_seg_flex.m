function res = check_same_seg_flex(seg1,seg2,num)
%check whether two segmentations are the same without assuming 
%that that when merging, the segment containing the smallest index will be 
%the one remained 
%input:
% seg1,seg2: 1 by num cell array representing two different segmentations,
% each cell represents a connected subsets
%output:
% True if the two segmentations are identical 


keySet = [];
for i = 1:num
    if isempty(seg1{i})
        continue
    else
        keySet = [keySet, strjoin(string(sort(seg1{i}))," ")];
    end       
end
M = containers.Map(keySet,zeros(1,length(keySet)));
res = true;
for i = 1:num
    if isempty(seg2{i})
        continue
    else
        currstr = strjoin(string(sort(seg2{i}))," ");
        if ~isKey(M,currstr)
            res = false;
            break
        else
            M(currstr) = 1;
        end
    end       
end

for i = 1:num
    if isempty(seg2{i})
        continue
    else
        currstr = strjoin(string(sort(seg2{i}))," ");
        if ~isKey(M,currstr)
            res = false;
            break
        else
            M(currstr) = 1;
        end
    end       
end
v = M.values;
if any([v{:}]==0)
    res = false;
end


