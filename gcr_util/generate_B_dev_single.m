function Bdev = generate_B_dev_single(t,M)
%generate the 4Mp+2 by 2 derivative of B(t) matrix 
a = div2(M-1);
for k = 1:M
    Bdev(k,:) = [-2*pi*(k-a-1)*sin(2*pi*t*(k-a-1)),2*pi*(k-a-1)*cos(2*pi*t*(k-a-1))];
end
for k = (M+1):2*M
    Bdev(k,:) = [-cos(2*pi*t*(k-a-M-1))*2*pi*(k-a-M-1),-sin(2*pi*t*(k-a-M-1))*2*pi*(k-a-M-1)];
end
Bdev(M+a+2,:)=[];

end

function x = div2(y)
    if y>0
        x = floor(y/2);
    else
        x = ceil(y/2);
    end
end
     