function [] = plot_segmentation_pause(DT, selected, cx, cy, colors,time,contour,frame,pt_size)
%plot the segmentation, pause between each one and display the segment's id
%in the middle of the region 
%frame: restricted the plot to certain region, default by [0 1 0 1]
GRAY = [0.6 0.6 0.6];
%figure
%triplot(DT, 'Color', GRAY)
hold on
% the final result
index = 0;
reg_in = 0;
reg_in_set=[];
for i = reshape(contour.regionId,1,[])

    if ~isempty(selected{i})
        index = index+1;
        X_pt = cx(selected{i});
        Y_pt = cy(selected{i});
        disp(length(X_pt))
                disp(length(Y_pt))

      %place the region's Id in the frame close to the center of the seg
        in_frame = X_pt>=frame(1) & X_pt<=frame(2) &Y_pt>=frame(3) & Y_pt<=frame(4);
        X_pt = X_pt(in_frame);
        Y_pt = Y_pt(in_frame);
        if isempty(X_pt) || isempty(Y_pt)
            continue
        end
        scatter(X_pt,Y_pt, pt_size,  colors(index, :), '.')

        reg_in = reg_in+1;
        reg_in_set(reg_in) = i;
        %         dist = (X_pt-mean(X_pt)).^2 +(Y_pt-mean(Y_pt)).^2; 
%         [~,idx_dist] = sort(dist);
        
        pause(time)

    end
        axis equal

    axis(frame)
end

for i = reshape(contour.regionId,1,[])
    currseg = contour.contourV{i};
        plot([currseg(:,1)'; currseg(:,2)'], [currseg(:,3)'; currseg(:,4)'],'k');
  
end
for i = reshape(contour.regionId,1,[])
        X_pt = cx(selected{i});
        Y_pt = cy(selected{i});
        in_frame = X_pt>=frame(1) & X_pt<=frame(2) &Y_pt>=frame(3) & Y_pt<=frame(4);
        X_pt = X_pt(in_frame);
        Y_pt = Y_pt(in_frame);
        if isempty(X_pt) || isempty(Y_pt)
            continue
        end
        min_idx = min_sum_dist(X_pt,Y_pt);
        text(X_pt(min_idx),Y_pt(min_idx),num2str(i),'FontSize',5,'FontWeight','bold');
    
end
axis(frame)
[leg,icons] = legend(string(reg_in_set));
 icons = findobj(icons,'Type','patch');
 icons = findobj(icons,'Marker','none','-xor');
 set(icons,'MarkerSize',20);

leg.FontSize=8;
set(leg,'AutoUpdate','off')
title(strcat('the number of segments in frame is ',num2str(reg_in)))

    %find the point that has smallest summed distance w.r.t. the rest 
    function min_idx = min_sum_dist(X_pt,Y_pt)
        for idx = 1:length(X_pt)
            abs_dist(idx) = sum(abs(X_pt-X_pt(idx))+abs(Y_pt-Y_pt(idx)));
        end
        [~,min_idx] = min(abs_dist);
    end


end
