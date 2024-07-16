function [x_sample,y_sample] = sample_curve(x,y,K,x_start,clockwise)
%Sample K point from a closed curve based on equal arclength.
%Input:
% x: X coordinates of the closed curve
% y: Y coordinates of the closed curve
% K: number of sample points 
% c: the x coordinate of the first points
% clockwise: logic argument, if true then sample in a clockwise manner 
%
%Output: 
% x_sample: X coordinates of the K sampled points along the closed curve 
% y_sample: Y coordinates of the K sampled points along the closed curve 
%


x_sample = zeros(K,1);
y_sample = zeros(K,1);
xlimit = [x_start x_start];
ylimit = [0 1];

if clockwise==true
    if ispolycw(x, y) 
    else 
        x = flip(x);
        y = flip(y);
    end
else 
    if ispolycw(x, y)
        x = flip(x);
        y = flip(y);
    else 

    end

end

%adjust x,y according to the starting point
[xi,yi,ii] = polyxpoly(x,y,xlimit,ylimit);
[y_start,i] = max(yi);
x = [x_start;x(ii(i)+1:length(x));x(1:ii(i))];
y = [y_start;y(ii(i)+1:length(y));y(1:ii(i))];
n = length(x);
arcL = zeros(1,n);
x = [x;x(1)];
y = [y;y(1)];





for i = 1:n
    arcL(i) = sqrt((x(i+1)-x(i))^2+(y(i+1)-y(i))^2);
end

arcL = [0, arcL];
cumL = cumsum(arcL);
diffL = cumL(n+1)/K;



currArc= [1,2];
x_sample(1) = x(1);
y_sample(1) = y(1);

for k = 1:(K-1)
    while k*diffL >= cumL(currArc(2)) 
        currArc = currArc+1;
    end
    L = k*diffL-cumL(currArc(1));
    x_sample(k+1) = x(currArc(1))+L/arcL(currArc(2))*(x(currArc(2))-x(currArc(1)));
    y_sample(k+1) = y(currArc(1))+L/arcL(currArc(2))*(y(currArc(2))-y(currArc(1)));  
end

end

