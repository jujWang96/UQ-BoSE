cd('/Users/joycewang/src/gsrg/plot_util')
addpath(genpath('~/src/gsrg/FD_util'))

close all

% plot an illustration figure for concentric shape 
unit_a = 0.03;
radius = [(unit_a/pi)^(1/2), (4*unit_a/pi)^(1/2),(10*unit_a/pi)^(1/2)];
figure;
plot_circles([0.5,0.52; 0.48,0.5;0.55,0.5], radius,1,'k')
axis equal
box on
axis([0,1,0,1])

set(gcf,'position',[50,50,300,300])
xticks(linspace(0,1,6))
yticks(linspace(0,1,6))
saveas(gcf, 'concentric2_seg.png');

% plot an illustration figure for the nonconcentric and intersecting curves
unit_a = 0.03;
radius = [(unit_a/pi)^(1/2), (4*unit_a/pi)^(1/2),(20*unit_a/pi)^(1/2)];
figure;
plot_circles([0.5,0.52; 0.48,0.5;0.6,0.5], radius,1,'k')
axis equal
box on
axis([0,1,0,1])

set(gcf,'position',[50,50,300,300])
xticks(linspace(0,1,6))
yticks(linspace(0,1,6))
saveas(gcf, 'nonconcentric2_seg.png');

% plot an illustration figure for the nonconcentric and intersecting curves
unit_a = 0.03;
radius = [(unit_a/pi)^(1/2), (4*unit_a/pi)^(1/2),(10*unit_a/pi)^(1/2)];
figure;
plot_circles([0.57772050238,0.5; 0.48,0.5;0.55,0.5], radius,1,'k')
axis equal
box on
axis([0,1,0,1])

set(gcf,'position',[50,50,300,300])
xticks(linspace(0,1,6))
yticks(linspace(0,1,6))
saveas(gcf, 'concentric3_seg.png');

% plot an illustration figure for correcting nonconcentric and intersecting curves
unit_a = 0.03;
radius = [(unit_a/pi)^(1/2), (4*unit_a/pi)^(1/2),(20*unit_a/pi)^(1/2)];
figure;
plot_circles([0.5,0.52; 0.48,0.5;0.6,0.5], radius,1,'k')
axis equal
box on
axis([0,1.21,0,1.21])

set(gcf,'position',[50,50,300,300])
xticks(linspace(0,1,6))
yticks(linspace(0,1,6))
saveas(gcf, 'correction2_1_seg.png');

% plot an illustration figure for correcting nonconcentric and intersecting curves
unit_a = 0.03;
radius = [(unit_a/pi)^(1/2), (4*unit_a/pi)^(1/2)];
figure;
plot_circles([0.5,0.52; 0.48,0.5], radius,1,'k')
axis equal
box on
axis([0,1,0,1])

set(gcf,'position',[50,50,300,300])
xticks(linspace(0,1,6))
yticks(linspace(0,1,6))
saveas(gcf, 'correction2_2_seg.png');


% plot an illustration figure for the concentric and non-intersecting curves
unit_a = 0.03;
radius = [(unit_a/pi)^(1/2),(4*unit_a/pi)^(1/2)];
figure;
plot_circles([0.3,0.3;0.6,0.6], radius,1,'k')
axis equal
box on
axis([0,1,0,1])

set(gcf,'position',[50,50,300,300])
xticks(linspace(0,1,6))
yticks(linspace(0,1,6))
saveas(gcf, 'nonconcentric3_seg.png');

% plot an illustration figure for the nonconcentric and intersecting curves
unit_a = 0.03;
radius = [(4*unit_a/pi)^(1/2),(4*unit_a/pi)^(1/2)];
figure;
plot_circles([0.45,0.5;0.55,0.5], radius,1,'k')
axis equal
box on
axis([0,1,0,1])

set(gcf,'position',[50,50,300,300])
xticks(linspace(0,1,6))
yticks(linspace(0,1,6))
saveas(gcf, 'nonconcentric1_seg.png');


[x_left, y_left] = generateCircleCoordinates((4*unit_a/pi)^(1/2), 0.45, 0.5, 100);
[x_right, y_right] = generateCircleCoordinates((4*unit_a/pi)^(1/2), 0.55, 0.5, 100);

pgon_left = polyshape(x_left, y_left);
pgon_right = polyshape(x_right,y_right);
[x_int, y_int] = polybool('intersection', x_left, y_left, x_right, y_right);
[x_diff_l, y_diff_l] = polybool('subtraction', x_left, y_left, x_right, y_right);
[x_diff_r, y_diff_r] = polybool('subtraction', x_right, y_right, x_left, y_left);

pgon_int = polyshape(x_int,y_int);
pgon_diff_l = polyshape(x_diff_l,y_diff_l);
pgon_diff_r = polyshape(x_diff_r,y_diff_r);

figure;
plot(pgon_diff_l)
hold on 
plot(pgon_diff_r)
plot(pgon_int)

axis equal
box on
axis([0,1,0,1])
set(gcf,'position',[50,50,300,300])
saveas(gcf, 'nonconcentric1_polygon.png');


%extract the boundary 
[x_boundary_l, y_boundary_l] = poly2cw(x_diff_l, y_diff_l);
[x_boundary_r, y_boundary_r] = poly2cw(x_diff_r, y_diff_r);
[x_boundary_int, y_boundary_int] = poly2cw(x_int, y_int);

figure 
lineWidth = 1;
plot (x_boundary_l-0.03,y_boundary_l, 'LineWidth', lineWidth)
hold on 
plot (x_boundary_r+0.03,y_boundary_r, 'LineWidth', lineWidth)
plot (x_boundary_int,y_boundary_int, 'LineWidth', lineWidth)
axis equal
box on
axis([0,1,0,1])
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gcf,'position',[50,50,300,300])
saveas(gcf, 'nonconcentric1_extract_curve.png');

coef_l = FD(x_boundary_l, y_boundary_l);
coef_r = FD(x_boundary_r, y_boundary_r);
coef_int = FD(x_boundary_int, y_boundary_int);
M = 15;
[x_rboundary_l, y_rboundary_l] = iFD(reshape(coef_l, [], 1) ,M);
[x_rboundary_r, y_rboundary_r] = iFD(reshape(coef_r, [], 1),M);
[x_rboundary_int, y_rboundary_int] = iFD(reshape(coef_int, [], 1),M+1);
figure
plot (x_rboundary_l,y_rboundary_l, 'LineWidth', lineWidth)
hold on 
plot (x_rboundary_r,y_rboundary_r, 'LineWidth', lineWidth)
plot (x_rboundary_int,y_rboundary_int, 'LineWidth', lineWidth)
axis equal
box on
axis([0,1,0,1])
set(gcf,'position',[50,50,300,300])
hold off
saveas(gcf, 'nonconcentric1_param_curve.png');

% 

rpgon_l = polyshape(x_rboundary_l, y_rboundary_l);
rpgon_r = polyshape(x_rboundary_r, y_rboundary_r);
rpgon_int = polyshape(x_rboundary_int, y_rboundary_int);

[distinctRegionsX, distinctRegionsY] = partitionPolygons({x_rboundary_l,x_rboundary_r,x_rboundary_int}, {y_rboundary_l,y_rboundary_r,y_rboundary_int});
figure
hold on 
for i = 1:length(distinctRegionsX)
    plot(polyshape(distinctRegionsX{i}, distinctRegionsY{i}))
end
axis equal
box on
axis([0,1,0,1])
set(gcf,'position',[50,50,300,300])
hold off
saveas(gcf, 'nonconcentric1_param_polygon.png');

figure 
lineWidth = 1;
plot (x_boundary_l,y_boundary_l, 'LineWidth', lineWidth, 'Color','k')
hold on 
plot (x_boundary_int,y_boundary_int, 'LineWidth', lineWidth, 'Color','k')
axis equal
box on
axis([0,1,0,1])
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gcf,'position',[50,50,300,300])
saveas(gcf, 'nonconcentric4_seg.png');


figure 
lineWidth = 1;
pgon_confine = polyshape(x_boundary_l,y_boundary_l);
plot (pgon_confine)
hold on 
plot (x_boundary_int,y_boundary_int, 'LineWidth', lineWidth, 'Color','r')
axis equal
box on
axis([0,1,0,1])
% set(gca, 'XTick', []);
% set(gca, 'YTick', []);
set(gcf,'position',[50,50,300,300])
saveas(gcf, 'violates_enclosure','epsc');

unit_a = 0.03;
radius = [(4*unit_a/pi)^(1/2)];
figure;
[x_left, y_left] = generateCircleCoordinates((4*unit_a/pi)^(1/2), 0.45, 0.5, 100);

plot(x_left, y_left,'LineWidth', lineWidth)
hold on 
plot (x_boundary_int-0.03,y_boundary_int, 'LineWidth', lineWidth)
axis equal
box on
axis([0,1,0,1])
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gcf,'position',[50,50,300,300])
saveas(gcf, 'nonconcentric4_extract_curve.png');

%extract the boundary 
[x_boundary_l, y_boundary_l] = poly2ccw(x_left, y_left);
[x_boundary_int, y_boundary_int] = poly2ccw(x_int, y_int);
coef_l = FD(x_boundary_l, y_boundary_l);
coef_int = FD(x_boundary_int, y_boundary_int);
M = 15;
[x_rboundary_l, y_rboundary_l] = iFD(reshape(coef_l, [], 1) ,M);
[x_rboundary_int, y_rboundary_int] = iFD(reshape(coef_int, [], 1),M+1);
figure
plot (x_rboundary_l,y_rboundary_l, 'LineWidth', lineWidth)
hold on 
plot (x_rboundary_int,y_rboundary_int, 'LineWidth', lineWidth)
axis equal
box on
axis([0,1,0,1])
set(gcf,'position',[50,50,300,300])
hold off
saveas(gcf, 'nonconcentric4_param_curve.png');

M = 2;
[x_rboundary_int, y_rboundary_int] = iFD(reshape(coef_int, [], 1),M);
figure
plot (pgon_confine)
hold on 
plot (x_rboundary_int,y_rboundary_int, 'LineWidth', lineWidth)
axis equal
box on
axis([0,1,0,1])
set(gcf,'position',[50,50,300,300])
hold off
saveas(gcf, 'nonconcentric4_model_selection_2fd','epsc');


[distinctRegionsX, distinctRegionsY] = partitionPolygons({x_rboundary_l,x_rboundary_int}, {y_rboundary_l,y_rboundary_int});
figure
hold on
for i = 1:length(distinctRegionsX)

    plot(polyshape(distinctRegionsX{i}, distinctRegionsY{i}))
end
axis equal
box on
axis([0.3,0.7,0.3,.7])
set(gcf,'position',[50,50,300,300])
hold off
saveas(gcf, 'nonconcentric4_param_polygon.png');



% plot an illustration figure for correcting nonconcentric and intersecting curves
unit_a = 0.03;
radius = [(unit_a/pi)^(1/2), (4*unit_a/pi)^(1/2)];
figure;
plot_circles([0.5,0.52; 0.48,0.5], radius,1,'k')
hold on 
axis equal
box on
axis([0,1.21,0,1.21])
x_line = linspace(0.5,1,100);
y_line = x_line*0.2+0.5;
plot(x_line,y_line,'LineWidth', lineWidth,'Color','k')
