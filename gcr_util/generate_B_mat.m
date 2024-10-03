function Bt = generate_B_mat(t,Mp)
%generate the g by 2 B(t) matrix 
for k = 1:(2*Mp+1)
    Bt(k,:) = [cos(2*pi*t*(k-Mp-1)),sin(2*pi*t*(k-Mp-1))];
end
for k = (2*Mp+2):(4*Mp+2)
    Bt(k,:) = [-sin(2*pi*t*(k-(2*Mp+1)-Mp-1)),cos(2*pi*t*(k-(2*Mp+1)-Mp-1))];
end
Bt(3*Mp+3,:)=[];
end


