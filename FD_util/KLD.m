function d =  KLD(X,coef0,coef1,p0in,p0out,p1in,p1out,M0,M1)
%compute the information lost if fitting X to model 1 than fitting to model
%0 i.e., D_{KL}() 

if nargin==7
    M0 = length(coef0);
    M1 = length(coef1);
end
% n0=length(X);
 cx = X(:,1);
 cy = X(:,2);
 [xv0,yv0] = iFD(coef0,M0);
%  A0 = polyarea(xv0,yv0);
 in0 = inpolygon(cx,cy,xv0,yv0);
[xv1,yv1] = iFD(coef1,M1);
 
 in1 = inpolygon(cx,cy,xv1,yv1);
% A1 = polyarea(xv1,yv1);
% l0in = sum(in0)/A0/n0;
% l0out = sum(~in0)/(1-A0)/n0;
% p1in = sum(in1)/A1/n0;
% p1out = sum(~in1)/(1-A1)/n0;
% if nargin==5
%     pin = p1in;
%     pout = p1out;
% end

d=0;
for i = 1:length(X)
    if in0(i)
        if in1(i)
            d = d+p0in*(log(p0in)-log(p1in));
        else
            d = d+p0in*(log(p0in)-log(p1out));
        end
        
    else
        if in1(i)
            d = d+p0out*(log(p0out)-log(p1in));
        else
            d = d+p0out*(log(p0out)-log(p1out));
        end
    end
end

% for i = 1:n0
%     if in0(i)
%         if in1(i)
%             d = d+l0in*(log(l0in)-log(l1in));
%         else
%             d = d+l0in*(log(l0in)-log(l1out));
%         end
%         
%     else
%         if in1(i)
%             d = d+l0out*(log(l0in)-log(l1in));
%         else
%             d = d+l0out*(log(l0in)-log(l1out));
%         end
%     end
% end
end