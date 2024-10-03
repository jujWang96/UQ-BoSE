function p = calc_lc_bd_single(L,M,C,beta)
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
 for l=1:L
     k0=k0+((calc_xi((l-0.5)/L,M,C)/2/pi)^(1/2))/L;
 end
 disp(k0)
 p = 1-(1+k0*beta)*exp(-beta^2/2);
%p =  calc_xi(0.9,Mp,C);
function xi  = calc_xi(t,M,C)
     B = generate_B_mat_single(t,M);
     dB = generate_B_dev_single(t,M);
    %B = generate_B_mat2(t)';
    %dB = generate_B_dev2(t)';
    Chf = sqrtm(C);
    A = (Chf*B)*B'*Chf;
    dA = Chf*(dB*B'+B*dB')*Chf;
    Ainv = pinv(A);
    n=length(A);
    Cf = B'*C*B;
    %disp(Chf*dB*pinv(Cf)*dB'*Chf)

    %xi = trace((eye(n)-A*pinv(A'*A)*A')*dA*Ainv*Ainv*dA*(eye(n)-A*pinv(A'*A)*A'));
    xi = trace((eye(n)-A*pinv(A'*A)*A')*Chf*dB*pinv(Cf)*dB'*Chf);
end
end