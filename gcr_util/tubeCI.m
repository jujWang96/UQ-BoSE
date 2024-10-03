function [CI,C] = tubeCI(origin_coef,boot_F_coef,L,Mp,p,r, polymask)
%construct the confidence region by moving ellipse along the curve 
%first extract the truncated FD
%p is the confidence level for each ellipse
%r is the resolution 
% displaying the confidence interval using imshow(flip(CI)) 
% M: number of FD **pair** used
% https://stats.stackexchange.com/questions/9898/how-to-plot-an-ellipse-from-eigenvalues-and-eigenvectors-in-r
if nargin==6
    polymask = true;
end
[~,N] = size(boot_F_coef);
for col = 1:N
    boot_coef_inv = make_startpt_inv2(boot_F_coef(:,col));

    coef_truncate_inv(:,col) = boot_coef_inv([end-Mp+1:end,1:Mp+1]);
end
%para = transpose([real(coef_truncate);imag(coef_truncate)]);
para = transpose([real(coef_truncate_inv);imag(coef_truncate_inv)]);
para(:,3*Mp+3)=[];
origin_coef_inv = make_startpt_inv2(origin_coef);
origin_coef_inv = origin_coef_inv([end-Mp+1:end,1:Mp+1]);
mu = [real(origin_coef_inv);imag(origin_coef_inv)];
mu(3*Mp+3)=[];
%C is a 4M_p+2 squared real matrix 
C = cov(para);


t = linspace(0, 2 * pi);
s=chi2inv(p, 2);
count = zeros(r,r);
if polymask
    %construct the CI using polymask
    for k=1:L
        Bt = generate_B_mat(k/L,Mp);
        C_f = transpose(Bt)*C*Bt;
        Bdev = generate_B_dev(k/L,Mp);
        [V, D] = eig(C_f * s);
    
        %tagt = transpose(Bdev)*mu;
        a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
        ct(k,:) = transpose(transpose(Bt)*mu);
    
        %plot(a(1, :) + ct(k,1), a(2, :) + ct(k,2),'b');
        count = count+double(poly2mask((a(1, :) + ct(k,1))*r,(a(2, :) + ct(k,2))*r,r,r));
    end
    CI = zeros(r,r);
    CI(count>0)=1;
else
    %contruct the CI using polyshape
    CI = polyshape();
    for k=1:L
        Bt = generate_B_mat(k/L,Mp);
        C_f = transpose(Bt)*C*Bt;
        Bdev = generate_B_dev(k/L,Mp);
        [V, D] = eig(C_f * s);
    
        %tagt = transpose(Bdev)*mu;
        a = (V * sqrt(D)) * [cos(t(:))'; sin(t(:))'];
        ct(k,:) = transpose(transpose(Bt)*mu);
    
        %plot(a(1, :) + ct(k,1), a(2, :) + ct(k,2),'b');
        CI = simplify(union(CI,polyshape((a(1, :) + ct(k,1)),(a(2, :) + ct(k,2)))));
    end
end
%-------------------------------------------------------------
% %generate the matrix B of 
%  function Bt = generate_B_mat(t,Mp)
% %generate the 4Mp+2 by 2 B(t) matrix 
% for k = 1:2*Mp+1
%     Bt(k,:) = [cos(2*pi*t*(k-Mp-1)),sin(2*pi*t*(k-Mp-1))];
% end
% for k = 2*Mp+2:4*Mp+2
% %     if k ==2*Mp+3
% %         Bt(k,:) = [0,0];
% %     else
%         Bt(k,:) = [-sin(2*pi*t*(k-(2*Mp+1)-Mp-1)),cos(2*pi*t*(k-(2*Mp+1)-Mp-1))];
% %     end
% end
% end
%---------------------------------------------------------------
% getpoints - Find the two tangent points to the ellipse with the given (dy/dx) slope, given a 2x2
%   covariance matrix, C_f and center mu, with e1 the closer one to
%   centroid
%
% function [e1,e2] = get_tangent_points(slope,C_f,mu,p,centrd)
% s=chi2inv(p, 2);
% [V, D] = eig(C_f * s);
% % [l, ind] = sort(l,'descend');
% % V = V(:, ind);
% l = diag(D);
% theta = atan(-sqrt(l(2))/sqrt(l(2))/slope);
% %plot(sqrt(l(1))*cos(rot)*cos(t)-sqrt(l(2))*sin(t)*sin(rot),sqrt(l(1))*sin(rot)*cos(t)+sqrt(l(2))*sin(t)*cos(rot));
% %error_ellipse(C_f);
% % axis equal
%     e1 = [sqrt(l(1))*cos(theta)*cos(0)-sqrt(l(2))*sin(0)*sin(theta),sqrt(l(1))*sin(theta)*cos(0)+sqrt(l(2))*sin(0)*cos(theta)]+mu;
%     e2 = [sqrt(l(1))*cos(theta)*cos(pi)-sqrt(l(2))*sin(pi)*sin(theta),sqrt(l(1))*sin(theta)*cos(pi)+sqrt(l(2))*sin(pi)*cos(theta)]+mu;
% 
% if norm(centrd-e1)>norm(centrd-e2)
%     [e2, e1] = deal(e1,e2);
% end
% % hold on 
% % scatter(e1(1),e1(2));
% % hold on 
% % scatter(e2(1),e2(2));
% end
% 
% %---------------------------------------------------------------
% % getpoints - Find the two end points on the major exis on the ellipse, given a 2x2
% %   covariance matrix, C_f and center mu, with e1 the closer one to
% %   centroid
% %
% function [e1,e2] = get_points(C_f,mu,p,centrd)
% s=chi2inv(p, 2);
% [V, D] = eig(C_f * s);
% % [l, ind] = sort(l,'descend');
% % V = V(:, ind);
% l = diag(D);
% %find the rotation of ellipse
% if C_f(1,2)==0 && C_f(1,1)>=C_f(2,2)
%     rot = 0;
% elseif C_f(1,2)==0 && C_f(1,1)<C_f(2,2)
%     rot = pi/2;
% else
%     rot = atan2(l(1)-s*C_f(1,1),s*C_f(1,2));
% end
% %plot(sqrt(l(1))*cos(rot)*cos(t)-sqrt(l(2))*sin(t)*sin(rot),sqrt(l(1))*sin(rot)*cos(t)+sqrt(l(2))*sin(t)*cos(rot));
% %error_ellipse(C_f);
% % axis equal
% clear e1 e1
% if l(1)>=l(2)
%     e1 = [sqrt(l(1))*cos(rot)*cos(0)-sqrt(l(2))*sin(0)*sin(rot),sqrt(l(1))*sin(rot)*cos(0)+sqrt(l(2))*sin(0)*cos(rot)]+mu;
%     e2 = [sqrt(l(1))*cos(rot)*cos(pi)-sqrt(l(2))*sin(pi)*sin(rot),sqrt(l(1))*sin(rot)*cos(pi)+sqrt(l(2))*sin(pi)*cos(rot)]+mu;
% else
%     e1 = [sqrt(l(1))*cos(rot)*cos(pi/2)-sqrt(l(2))*sin(pi/2)*sin(rot),sqrt(l(1))*sin(rot)*cos(pi/2)+sqrt(l(2))*sin(pi/2)*cos(rot)]+mu;
%     e2 = [sqrt(l(1))*cos(rot)*cos(3*pi/2)-sqrt(l(2))*sin(3*pi/2)*sin(rot),sqrt(l(1))*sin(rot)*cos(3*pi/2)+sqrt(l(2))*sin(3*pi/2)*cos(rot)]+mu;
% 
% end
% if norm(centrd-e1)>norm(centrd-e2)
%     [e2, e1] = deal(e1,e2);
% end
% % hold on 
% % scatter(e1(1),e1(2));
% % hold on 
% % scatter(e2(1),e2(2));
% end

end