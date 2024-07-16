function [] = gSRG_reproduce(X,glb_contour,glb_region_intensity,glb_selected,seg_contour,seg_region_intensity,seg_selected,frame) 
%
%
%Use SRG-graph to segment the points, merge them uniform randomly and
%return the contour of the segmentation. 
%Input: 
% X: all of the points
% draw_s: logic; if true then draw the segmentation
% draw_c: logic; if true then draw the contour
% Color: the color of contour
% N_r: the true number of region 
% drop: logic, if true then drop this simulation when the segmentation 
% not equal to true numbe of region
% restrict: if the segmentation region number is less than restrict number, 
% return segmentation of N_r regions 
% rand_num: number of possible merge to choose each time
% rand_itr: the iteration that uses random merge
% m: penalty coefficient
%
%Output:
% C: cell array, each element contains x1,x2,y1,y2 representing the 
% contour of of one region 

subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.05 0.02], [0.05 0.02]);


% init comp
[cx, cy, n, DT, E, cell_log_intensity, cell_area] = init_comp(X, [0 1], [0 1], ones(size(X, 1), 1));
num = length(glb_region_intensity);

subplot(2,2,2)
plot_intensity( cx,cy, glb_selected, glb_contour,glb_region_intensity,0.3)
axis(frame)
[leg,icons] = legend(string(glb_contour.regionId));
 icons = findobj(icons,'Type','patch');
 icons = findobj(icons,'Marker','none','-xor');
 set(icons,'MarkerSize',15);
leg.FontSize = 8;
axis(frame)
subplot(2,2,1)
plot_segmentation_pause(DT, glb_selected, cx, cy, lines(num),0,glb_contour,frame,0.3)


subplot(2,2,3)
plot_segmentation_pause(DT, seg_selected, cx, cy, lines(num),0,seg_contour,frame,0.3)

set(gcf, 'Position', [50 50 600 600]); %

end


