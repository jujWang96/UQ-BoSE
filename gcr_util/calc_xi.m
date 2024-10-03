function xi  = calc_xi(t,Mp,C)
     B = generate_B_mat(t,Mp);
     dB = generate_B_dev(t,Mp);
    %B = generate_B_mat2(t)';
    %dB = generate_B_dev2(t)';
    Chf = sqrtm(C);
    A = (Chf*B)*B'*Chf;
    %dA = Chf*(dB*B'+B*dB')*Chf;
    %Ainv = pinv(A);
    n=length(A);
    Cf = B'*C*B;
    %disp(Chf*dB*pinv(Cf)*dB'*Chf)

    %xi = trace((eye(n)-A*pinv(A'*A)*A')*dA*Ainv*Ainv*dA*(eye(n)-A*pinv(A'*A)*A'));
    xi = (trace((eye(n)-A*pinv(A'*A)*A')*Chf*dB*pinv(Cf)*dB'*Chf)/2/pi)^(1/2);
end
