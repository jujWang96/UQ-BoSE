function coef = get_coef(X)
%get the Fourier coefficients of the contour; 
K=300;
[x1,x2,y1,y2] = segment_graph_contour(X,false,false,'r',2,true,inf);
if isempty(x1) 
    coef = [];
    return;
end

[x,y] = get_curve(x1,x2,y1,y2,true);
[centx,centy] = centroid(polyshape(x,y));
[x_sample,y_sample] = sample_curve(x,y,K,centx, false);
coef=FD(x_sample,y_sample);
end
