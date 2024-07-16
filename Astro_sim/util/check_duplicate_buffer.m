function res = check_duplicate_buffer(buffer,bw)
%check whether there are duplicate segmentation stored in a buffer 
%return false if there is no


segstrSet = {};
res = false;

for jj = 1:bw
    keySet = [];
    seg = buffer{jj};
    num = length(buffer{jj});

    for ii = 1:num
        if isempty(seg{ii})
            continue
        else
            keySet = [keySet, strjoin(string(sort(seg{ii}))," ")];
        end       
    end
    segstr = strjoin(sort(keySet),",");
    %disp(segstr)
    if ismember(segstr,segstrSet)
        res = true;

        break
    else
        segstrSet = [segstrSet,segstr];
    end
end
