function p = calc_lc_bd_single_continuous(M,C,Beta)
%estimate the level-crossing probability 
%Input variables:
%
%L: number of samples along the curve
%Mp: number of paired FD to parametrize the curve (Mp = M//2)
%C: covariance matrix of parameters
%beta: level crossing probability at single point
%
%Output variable:
%
%p: The empirical level crossing probability  
k0=0;
Chf = sqrtm(C);

 k0 = integral(@(t) (calc_xi(t,M,C,Chf)/2/pi)^(1/2), 0,1,"ArrayValued",true);
 
 disp(k0)
 p = 1-(1+k0*Beta).*exp(-Beta.^2/2);
%p =  calc_xi(0.9,Mp,C);
function xi  = calc_xi(t,M,C,Chf)
     B = generate_B_mat_single(t,M);
     dB = generate_B_dev_single(t,M);
    %B = generate_B_mat2(t)';
    %dB = generate_B_dev2(t)';
    A = (Chf*B)*B'*Chf;
    
    n=size(A,1);
    Cf = B'*C*B;
    %disp(Chf*dB*pinv(Cf)*dB'*Chf)

    %xi = trace((eye(n)-A*pinv(A'*A)*A')*dA*Ainv*Ainv*dA*(eye(n)-A*pinv(A'*A)*A'));
    xi = trace((eye(n)-A*pinv(A'*A)*A')*Chf*dB*pinv(Cf)*dB'*Chf);
end
end