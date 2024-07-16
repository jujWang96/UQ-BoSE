function coef_res = rescale_coef(coef,x,y,dx,dy)
%convert fourier descriptor from [0,1]x[0,1] to a new scale proportially 
coef_res = coef;
for i = 1:size(coef,2)
    [sample_x,sample_y] = iFD(coef(:,i));
    rescale_x = sample_x*10*dx;
    rescale_y = sample_y*10*dy;
    coef_res(:,i) = FD(rescale_x,rescale_y);
end