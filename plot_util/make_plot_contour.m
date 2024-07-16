

function [] = make_plot_contour(X,contour,pt_size,line_size,target)
%target: target contour of length d represented by 2xd matrix
GREY = [0.6,0.6,0.6];

scatter(X(:,1),X(:,2),pt_size,  GREY, '.')
hold on
for c = contour.regionId
    contourv = contour.contourV{c};
    for d = 1:length(contourv)
        plot([contourv(d,1),contourv(d,2)],[contourv(d,3),contourv(d,4)],'b',LineWidth=line_size)
    end    
end
if nargin==5
    if size(target,2)==4
        for d = 1:length(target)
            
            plot([target(d,1),target(d,2)],[target(d,3),target(d,4)],'r',LineWidth=line_size)
            
        end
    else
        plot(target(:,1),target(:,2),'r',LineWidth=2)
    end
end
hold off
axis equal
axis([0 1 0 1])
set(gcf,'renderer','Painters')
box on 
axis equal
axis([0 1 0 1])