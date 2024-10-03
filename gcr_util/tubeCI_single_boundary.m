function [CI,C,boundaries] = tubeCI_single_boundary(origin_coef,boot_F_coef,L,M,p,r,polymask)
%construct the confidence region by moving ellipse along the curve 
%first extract the truncated FD of length M
%p is the confidence level for each ellipse
%r is the resolution 
% if polymask is true, CI will be a polymask scaled to r-by-r frame, 
%   display the confidence interval using imshow(flip(CI)) 
% if polymask if false, CI will be a polyshape in 1-by-1 frame,
%   display the confidence interval using plot(CI)
% boundary: the cells of boundaries 
%M: number of FD used
if nargin==6
    polymask = true;
end
[~,N] = size(boot_F_coef);
coef_truncate_inv = zeros(M,N);
for col = 1:N
    boot_coef_inv = make_startpt_inv2(boot_F_coef(:,col));
    coef_truncate_inv(:,col) = boot_coef_inv([end-div2(M-1)+1:end,1:div2(M)+1]);
end
%para = transpose([real(coef_truncate);imag(coef_truncate)]);
para = transpose([real(coef_truncate_inv);imag(coef_truncate_inv)]);
para(:,M+div2(M-1)+2)=[];
origin_coef_inv = make_startpt_inv2(origin_coef);
origin_coef_inv = origin_coef_inv([end-div2(M-1)+1:end,1:div2(M)+1]);
mu = [real(origin_coef_inv);imag(origin_coef_inv)];
mu(M+div2(M-1)+2)=[];
%C is a 2M-1 squared real matrix 
C = cov(para);


t = linspace(0, 2 * pi);
s=chi2inv(p, 2);
if polymask
    count = zeros(r,r);
    for k=1:L
        Bt = generate_B_mat_single(k/L,M);
        C_f = transpose(Bt)*C*Bt;
        [V, D] = eig(C_f * s);
        a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
        ct(k,:) = transpose(transpose(Bt)*mu);
        count = count+double(poly2mask((a(1, :) + ct(k,1))*r,(a(2, :) + ct(k,2))*r,r,r));
    end
    
    CI = zeros(r,r);
    CI(count>0)=1;
    
    boundaries_bm = bwboundaries(CI);
    %rescale the boundaries and switch x and y axis 
    boundaries = cell(length(boundaries_bm),1);
    for ii = 1:length(boundaries_bm)
        boundary_bm = boundaries_bm{ii};
        boundaries{ii} = [boundary_bm(:,2)/r,boundary_bm(:,1)/r];
    end
else
    warning('off')
    CI = polyshape();
    for k=1:L
        Bt = generate_B_mat_single(k/L,M);
        C_f = transpose(Bt)*C*Bt;
        [V, D] = eig(C_f * s);
        a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
        ct(k,:) = transpose(transpose(Bt)*mu);
        CI = simplify(union(CI,polyshape((a(1, :) + ct(k,1)),(a(2, :) + ct(k,2)))));
    end
    nbound = numboundaries(CI);
    boundaries = cell(nbound,1);
    for ii = 1: nbound
        boundaries{ii} = boundary(CI,ii);
    end
end
%--------------------------------------------------------------
%integer division function 
function x = div2(y)
    if y>0
        x = floor(y/2);
    else
        x = ceil(y/2);
    end
end
     
end