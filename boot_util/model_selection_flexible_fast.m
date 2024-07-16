function [candidate_aic,candidate_bic,pgon_confine] = model_selection_flexible_fast(X,target_x, target_y, seg_contour,candidates, plot_results)
%perform the OPTIMIZED model selection with parallel for FDs of single target curve assuming that the target
%contours will not intersect with the inner bound or outer bound. 
%**could deal with model than one disconnected boundries in each stage
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
pgon_confine = polyshape();
x_curves = {};
y_curves = {};
pgon_list = {};
outer_enclose_id = 0;
outer_match_idx = 0;
inner_enclose_id = 0;
inner_match_idx = 0;
n_sub = 0;
% identify the two segments that fully enclose the targer curve 
for objId = seg_contour.regionId
    [cxs,cys] = get_curve_flexible(seg_contour.contourV{objId},false);
    pgon_i = polyshape(cxs,cys);
    max_area_i = 0;
    pgon_list{end+1} =pgon_i;
    x_curves{end+1}= cxs;
    y_curves{end+1} = cys;
    is_match = false;
    matched_idx = 0;
    for i =1:length(cxs)
        max_area_i = max(max_area_i, polyarea(cxs{i},cys{i}));
        is_match_i = isequal(sort([cxs{i},cys{i}], 1), sort(pgon_target.Vertices, 1));
        if is_match_i
            matched_idx = i;
        end
        is_match = is_match || is_match_i;
    end
    if is_match && (round(max_area_i,3)> round(area_target,3))
        outer_enclose_id = objId;
        outer_match_idx = matched_idx;
        n_sub = n_sub + sum(isinterior(pgon_i,X(:,1),X(:,2)));
        pgon_confine = union(pgon_confine, pgon_i);
        pgon_confine = simplify(pgon_confine);
    end
    if is_match && (round(max_area_i,3) == round(area_target,3))
        inner_enclose_id = objId;
        inner_match_idx = matched_idx;
        n_sub = n_sub + sum(isinterior(pgon_i,X(:,1),X(:,2)));
        pgon_confine = union(pgon_confine, pgon_i);
        pgon_confine = simplify(pgon_confine);
    end
end

   
%     candidate_aic=0; candidate_bic=0;
%     return 
if inner_enclose_id==0
    disp("concentric condition is not satisfied for inner bound, look for a flexibile substitude")
    [bx,by, inside_pgon] = find_subsititude_inner(X,pgon_target, pgon_list);
%     figure
%     plot(subtract( pgon_target,inside_pgon))
    pgon_confine = union(pgon_confine, subtract( pgon_target,inside_pgon));
    pgon_confine = simplify(pgon_confine);
    x_curves{end+1}={bx;target_x};
    y_curves{end+1}={by;target_y};
    inner_enclose_id = length(x_curves);
    inner_match_idx = 2;
    n_sub = n_sub+num_target-sum(isinterior(inside_pgon,X(:,1),X(:,2)));
end

if outer_enclose_id==0
    disp("concentric condition is not satisfied for outer bound, look for a flexibile substitude")
    [bx,by,outside_pgon] = find_subsititude_outer(X,pgon_target, pgon_list);
%     figure
%     plot(subtract(outside_pgon, pgon_target))
    pgon_confine = union(pgon_confine, subtract(outside_pgon, pgon_target));
    pgon_confine = simplify(pgon_confine);
    x_curves{end+1}={bx;target_x};
    y_curves{end+1}={by;target_y};
    outer_enclose_id = length(x_curves);
    outer_match_idx = 2;
    n_sub = n_sub+sum(isinterior(outside_pgon,X(:,1),X(:,2)))-num_target;

end
disp(n_sub)
%perform model selection 
[centx,centy] = centroid(polyshape(target_x,target_y));
[x_sample,y_sample] = sample_curve(target_x,target_y,300,centx, false);
f_coefs = FD(x_sample,y_sample);

aic_vals = NaN(1, length(candidates));
bic_vals = NaN(1, length(candidates));
log_vals = NaN(1, length(candidates));

for i = 1:length(candidates)
    candidate = candidates(i);
    [invx,invy] = iFD(f_coefs,candidate);
    
    x_outer_enclose = x_curves{outer_enclose_id};
    x_outer_enclose{outer_match_idx} = invx;
    y_outer_enclose = y_curves{outer_enclose_id};
    y_outer_enclose{outer_match_idx} = invy;
    
    pgon_outer = polyshape(x_outer_enclose,y_outer_enclose);
    
    x_inner_enclose = x_curves{inner_enclose_id};
    x_inner_enclose{inner_match_idx} = invx;
    y_inner_enclose = y_curves{inner_enclose_id};
    y_inner_enclose{inner_match_idx} = invy;
    
    pgon_inner = polyshape(x_inner_enclose,y_inner_enclose);
    if pgon_outer.NumRegions>1 || pgon_inner.NumRegions>1
       disp("disgard for intersecting bound") 
       continue
    end
    figure
    plot(pgon_inner)
    hold on
    plot(pgon_outer)
    a_in = area(pgon_inner);
    n_in = sum(isinterior(pgon_inner,X(:,1),X(:,2)));
    a_out = area(pgon_outer);
    n_out = sum(isinterior(pgon_outer,X(:,1),X(:,2)));
    log_like = n_in*log(n_in/a_in)+n_out*log(n_out/a_out)-sum(log(1:n_sub))-n_sub;
    aic_vals(i) = -2*log_like+(candidate*2+1)*2;
    bic_vals(i) = -2*log_like+(candidate*2+1)*log(n_sub);
    log_vals(i) = log_like;
end
[min_aic,aic_idx] = min(aic_vals);
[min_bic,bic_idx] = min(bic_vals);
candidate_aic = candidates(aic_idx);
candidate_bic = candidates(bic_idx);

if plot_results
    figure
    line(candidates, aic_vals,'DisplayName', 'aic')
    line(candidates, bic_vals,'DisplayName', 'bic')
    line(candidates, -2*log_vals, 'DisplayName', 'log')
    xlabel('number of parameters');
    ylabel('log');
    legend('Location', 'Best'); % 'Best' places the legend in the best available position

end

end
function [bx,by,inside_pgon] = find_subsititude_inner(X,pgon_target, pgon_list)
    %merge segments inside target_x, target_y. 
    inside_pgon = polyshape();
    target_in = isinterior(pgon_target,X(:,1),X(:,2));
    target_in = find(target_in ~= 0);
    for i = 1:length(pgon_list)
        pgon_i = pgon_list{i};
        %check whether pgon_i is inside pgon_target
        %is_inside = all(inpolygon(pgon_i.Vertices(:, 1), pgon_i.Vertices(:, 2), pgon_target.Vertices(:, 1), pgon_target.Vertices(:, 2)));
        pgon_i_in =  isinterior(pgon_i,X(:,1),X(:,2));
        pgon_i_in = find(pgon_i_in ~= 0);
        is_inside = intersect(pgon_i_in, target_in);
         %empty intersection means the segment is outside target
        if isempty(is_inside)
            continue
        end

%         %check whether pgon_i is inside pgon_target by comparing the area
%         %of their intersection
%         intersection_polygon = intersect(pgon_i, pgon_target);
%         
%         if round(area(intersection_polygon),3) < round(area(pgon_i))
%             continue
%         end
%         %the 
% 
%         %check the intersection of vertex of small 
        intersection_vertices = intersect(pgon_target.Vertices, pgon_i.Vertices, 'rows');
        if ~isempty(intersection_vertices)
            continue
        end
        
        inside_pgon = union(inside_pgon, pgon_i);
        inside_pgon = simplify(inside_pgon);
        
    end
    figure
    plot(inside_pgon)
    [bx,by] = boundary(inside_pgon);
end

function [bx,by,outside_pgon_reverse] = find_subsititude_outer(X,pgon_target,pgon_list)
    %merge segments inside target_x, target_y. 
    outside_pgon = polyshape();
    target_in = isinterior(pgon_target,X(:,1),X(:,2));
    target_in = find(target_in ~= 0);
    for i = 1:length(pgon_list)
        pgon_i = pgon_list{i};
        pgon_i_in =  isinterior(pgon_i,X(:,1),X(:,2));
        pgon_i_in = find(pgon_i_in ~= 0);
        is_inside = intersect(pgon_i_in, target_in);
         %empty intersection means the segment is outside target
        if ~isempty(is_inside)
            continue
        end
        intersection_vertices = intersect(pgon_target.Vertices, pgon_i.Vertices, 'rows');
        if ~isempty(intersection_vertices)
            continue
        end
        outside_pgon = union(outside_pgon, pgon_i);
        outside_pgon = simplify(outside_pgon);
        
    end
    figure
    plot(outside_pgon)
    outside_pgon_reverse = holes(outside_pgon);
    [bx,by] = boundary(outside_pgon_reverse);
end