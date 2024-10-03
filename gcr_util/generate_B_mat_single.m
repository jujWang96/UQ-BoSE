function Bt = generate_B_mat_single(t,M)
%generate the g by 2 B(t) matrix 
a = div2(M-1);
for k = 1:M
    Bt(k,:) = [cos(2*pi*t*(k-a-1)),sin(2*pi*t*(k-a-1))];
end
for k = (M+1):2*M
    Bt(k,:) = [-sin(2*pi*t*(k-M-a-1)),cos(2*pi*t*(k-M-a-1))];
end
Bt(M+a+2,:)=[];
end


function x = div2(y)
    if y>0
        x = floor(y/2);
    else
        x = ceil(y/2);
    end
end
     