function [candidate_aic, candidate_bic, pgon_confine] = model_selection_flexible_robust(X,target_x, target_y, seg_contour,candidates, plot_results)
%perform the robust model selection for FDs of single target curve assuming that the target
%contours will not intersect with the inner bound or outer bound. 
%final version implementation
% the confined region is the union of all segments that share its boundary
% with the target curve (target_x, target_y)
%input:
%   X: 2d photons in fov
%   target_x: 1-by-d vector of x coordinates
%   target_y: 1-by-d vector of y coordinates
%   seg_contour: segment contour contains the line segments of each contour
%   of the segmentation of fov
%   candidates: the range of J

pgon_target = polyshape(target_x,target_y);
area_target = area(pgon_target);
num_target= sum(isinterior(pgon_target,X(:,1),X(:,2)));
pgon_list = {};
n_sub = 0;

% get polygon from each segment
for objId = seg_contour.regionId
    [cxs,cys] = get_curve_flexible(seg_contour.contourV{objId},false);
    pgon_i = polyshape(cxs,cys);
    pgon_list{end+1} = pgon_i;
end


target_in = isinterior(pgon_target,X(:,1),X(:,2));
target_in = find(target_in ~= 0);


[in_x_curves, in_y_curves,inner_match_idx] = find_confine_region(X,target_x,target_y, pgon_target,target_in,pgon_list, true );
[out_x_curves, out_y_curves, outer_match_idx] = find_confine_region(X,target_x,target_y, pgon_target,target_in,pgon_list, false );
% out_x_curves = {target_x,[0,1,1,0]'};out_y_curves = {target_y,[0,0,1,1]'};
% outer_match_idx = 1;
pgon_confine = simplify(union(polyshape(in_x_curves, in_y_curves),polyshape(out_x_curves, out_y_curves)));
n_sub = sum(isinterior(pgon_confine,X(:,1),X(:,2)));
[centx,centy] = centroid(polyshape(target_x,target_y));
[x_sample,y_sample] = sample_curve(target_x,target_y,300,centx, false);
f_coefs = FD(x_sample,y_sample);

aic_vals = NaN(1, length(candidates));
bic_vals = NaN(1, length(candidates));
exp_vals = NaN(1, length(candidates));

log_vals = NaN(1, length(candidates));

for i = 1:length(candidates)
    candidate = candidates(i);
    [invx,invy] = iFD(f_coefs,candidate);
    
    x_outer_enclose = out_x_curves;
    x_outer_enclose{outer_match_idx} = invx;
    y_outer_enclose = out_y_curves;
    y_outer_enclose{outer_match_idx} = invy;
    
    pgon_outer = polyshape(x_outer_enclose,y_outer_enclose);
    
    x_inner_enclose = in_x_curves;
    x_inner_enclose{inner_match_idx} = invx;
    y_inner_enclose = in_y_curves;
    y_inner_enclose{inner_match_idx} = invy;
    
    pgon_inner = polyshape(x_inner_enclose,y_inner_enclose);
    if pgon_outer.NumRegions>1 || pgon_inner.NumRegions>1
       %disp("disgard for intersecting bound") 
       continue
    end
%     figure
%     plot(pgon_inner)
%     hold on
%     plot(pgon_outer)
    a_in = area(pgon_inner);
    n_in = sum(isinterior(pgon_inner,X(:,1),X(:,2)));
    a_out = area(pgon_outer);
    n_out = sum(isinterior(pgon_outer,X(:,1),X(:,2)));
    log_like = n_in*log(n_in/a_in)+n_out*log(n_out/a_out)-sum(log(1:n_sub))-n_sub;
    aic_vals(i) = -2*log_like+(candidate*2+1)*2;
    bic_vals(i) = -2*log_like+(candidate*2+1)*log(n_sub);
    exp_vals(i) = -2*log_like+(candidate+2)*log(n_sub); %this is not penalizing enough, select too many when n increases

    log_vals(i) = log_like;
end

%writematrix([candidates;log_vals;bic_vals]','debug_robust_subset.txt','Delimiter','tab')
if all(isnan(log_vals))
    disp("Bug in model selection.")
end

[min_aic,aic_idx] = min(aic_vals);
[min_bic,bic_idx] = min(bic_vals);
candidate_aic = candidates(aic_idx);
candidate_bic = candidates(bic_idx);
if plot_results
    figure
    disp('number of photon inside the confined shape is')
    disp(n_sub)
    line(candidates, aic_vals,'DisplayName', 'aic','Color', 'red')
    line(candidates, bic_vals,'DisplayName', 'bic','Color', 'blue')
    line(candidates, -2*log_vals, 'DisplayName', 'log','Color', 'green')
    xlabel('number of parameters');
    ylabel('log');
    legend('Location', 'Best'); % 'Best' places the legend in the best available position

end

end
function [x_curves, y_curves,match_index] = find_confine_region(X,target_x,target_y,pgon_target,target_in,pgon_list, inside )
    pgon_confine = polyshape();
    x_curves = {};
    y_curves = {};
    for i = 1:length(pgon_list)
        pgon_i = pgon_list{i};
        %check whether pgon_i is inside pgon_target
        %is_inside = all(inpolygon(pgon_i.Vertices(:, 1), pgon_i.Vertices(:, 2), pgon_target.Vertices(:, 1), pgon_target.Vertices(:, 2)));
        pgon_i_in =  isinterior(pgon_i,X(:,1),X(:,2));
        pgon_i_in = find(pgon_i_in ~= 0);
        is_inside = intersect(pgon_i_in, target_in);
         %empty intersection means the segment is outside target
        
        if inside && isempty(is_inside)
            continue
        end
        if (~inside) && (~isempty(is_inside))
            continue
        end

       intersection_vertices = intersect(pgon_target.Vertices, pgon_i.Vertices, 'rows');
        if isempty(intersection_vertices)
            continue
        end
        
        pgon_confine = union(pgon_confine, pgon_i);
        pgon_confine = simplify(pgon_confine);

    end
%      figure
%      plot(pgon_confine)


    % in all the boundaries of pgon_confine 
    % find the index of boundary with which is the target curve
    num_b = numboundaries(pgon_confine);
    match_index = 0;
    sort_target_xy = sort([target_x,target_y], 1);
    for b = 1:num_b
        [bx,by] = boundary(pgon_confine,b);
        sort_bxy = sort([bx(1:(end-1)),by(1:(end-1))], 1);
        intersection = ismember(sort_target_xy, sort_bxy, 'rows');
        %is_match_b = isequal(sort(round([bx(1:(end-1)),by(1:(end-1))],3), 1), sort(round([target_x,target_y],3), 1));
        is_match_b = any(intersection);

        if is_match_b
            match_index = b;
        end
        x_curves{b} = bx;
        y_curves{b} = by;
    end
end