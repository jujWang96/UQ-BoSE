function same_set = check_same_buffer(buffer1,buffer2,bw)
%check whether buffer1 and buffer2 contain the same segmentation set
%return the number of segmentations appearing in both buffers


segstrSet1 = {};
for jj = 1:bw
    keySet = [];
    seg = buffer1{jj};
    num = length(buffer1{jj});

    for ii = 1:num
        if isempty(seg{ii})
            continue
        else
            keySet = [keySet, strjoin(string(sort(seg{ii}))," ")];
        end       
    end
    segstr = strjoin(sort(keySet),",");
    %disp(segstr)
    if ~ismember(segstr, segstrSet1)
        segstrSet1 = [segstrSet1,segstr];
    end

end


segstrSet2 = {};
same_set = 0;
for jj = 1:bw
    keySet = [];
    seg = buffer1{jj};
    num = length(buffer1{jj});

    for ii = 1:num
        if isempty(seg{ii})
            continue
        else
            keySet = [keySet, strjoin(string(sort(seg{ii}))," ")];
        end       
    end
    segstr = strjoin(sort(keySet),",");
    %disp(segstr)
    if ~ismember(segstr, segstrSet2)
        segstrSet2 = [segstrSet2,segstr];
    end

end
same_set = length(intersect(segstrSet1,segstrSet2))*2-length(segstrSet1)-length(segstrSet2);