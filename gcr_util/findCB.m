function CB = findCB(hist_count,a)
%given the histogram count, find the 1-a% confidence band corresponding to
%smallest area 
[r,c] = size(hist_count);
total = sum(hist_count(:));
CB = zeros(r,c);
M=max(hist_count(:));
[cx,cy]=find(hist_count==M);
volumn=0;
neighbour = containers.Map(max(hist_count(:)),[cx(1),cy(1)]);
adj = zeros(r,c);
loop=0;

while volumn<=total*(1-a)
    loop = loop+ 1;
    val = max(cell2mat(keys(neighbour)));
%     if val==0
%         break
%     end
    cord = neighbour(val);
    remove(neighbour, val);
    %volumn = volumn+val*length(cord);
    cx = cord(:,1);
    cy = cord(:,2);
    for i = 1:length(cx)
        if CB(cx(i),cy(i))==0
       
            adj(cx(i),cy(i))=1;
        end
        
    end
    %all pixels should have no gap in between, so give the pixels in gap
    %highest priority in neighbour map.

    for i = 1:length(cx)
        if CB(cx(i),cy(i))==1
            continue
        end
        CB(cx(i),cy(i))=1; 
        volumn=volumn+hist_count(cx(i),cy(i));
        %add/update the neighbours of the new added coordinates, skip if
        %the coordinates is already in CB
        if adj(cx(i)+1,cy(i))==0 || CB(cx(i)+2,cy(i))==1
            adj(cx(i)+1,cy(i))=1;
            if CB(cx(i)+2,cy(i))==1
                currval = hist_count(cx(i)+1,cy(i))+M;
            else
                currval = hist_count(cx(i)+1,cy(i));
            end
            %append the coordinate if key alreay exist, add new key if not 
            if isKey(neighbour,currval)
                curridx = neighbour(currval);
                curridx(end+1,:)=[cx(i)+1,cy(i)];

                neighbour(currval)=curridx;
            else
                neighbour(currval)= [cx(i)+1,cy(i)];
            end
        end
      
        if adj(cx(i)-1,cy(i))==0 ||CB(cx(i)-2,cy(i))==1
            adj(cx(i)-1,cy(i))=1;
            if CB(cx(i)-2,cy(i))==1
                currval = hist_count(cx(i)-1,cy(i))+M;
            else
                currval = hist_count(cx(i)-1,cy(i));
            end
            if isKey(neighbour,currval)
                curridx = neighbour(currval);
                curridx(end+1,:)=[cx(i)-1,cy(i)];
                neighbour(currval)=curridx;
            else
                 neighbour(currval)= [cx(i)-1,cy(i)];
            end
        end
        if adj(cx(i),cy(i)+1)==0 || CB(cx(i),cy(i)+2)==1
            adj(cx(i),cy(i)+1)=1;
            if CB(cx(i),cy(i)+2)==1
                currval = hist_count(cx(i),cy(i)+1)+M;
            else
                currval = hist_count(cx(i),cy(i)+1);
            end
            if isKey(neighbour,currval)
                curridx = neighbour(currval);
                curridx(end+1,:)=[cx(i),cy(i)+1];

                neighbour(currval)=curridx;
            else
                neighbour(currval)= [cx(i),cy(i)+1];
            end
        end
        if adj(cx(i),cy(i)-1)==0 || CB(cx(i),cy(i)-2)==1
            adj(cx(i),cy(i)-1)=1;
            if CB(cx(i),cy(i)-2)==1
                currval = hist_count(cx(i),cy(i)-1)+M;
            else
                currval = hist_count(cx(i),cy(i)-1);
            end
            if isKey(neighbour,currval)
                curridx = neighbour(currval);
                curridx(end+1,:)=[cx(i),cy(i)-1];

                neighbour(currval)=curridx;
            else
                 neighbour(currval)= [cx(i),cy(i)-1];
            end
        end
        
    end
%     V = values(neighbour);
%     for v = 1:length(V)
%         disp(V{v});
%     end
end
