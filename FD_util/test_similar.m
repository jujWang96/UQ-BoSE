function rej = test_similar(Xb,Xe,B)
%test whether the two observed data set come from the same distribution 
%output:
%   rej:true if they do not come from the same distribution 
if nargin==2
    B = 200;
end
LRratio = zeros(1,B);
X_bind = [Xb;Xe];
coef0 = get_coef(X_bind);
coefb = get_coef(Xb);
coefe = get_coef(Xe);
l=0;
rej = false;
for b = 1:B
    rng(b);  
    idx_b = randperm(length(Xb),floor(length(Xb)/2));
    idx_e = randperm(length(Xe),ceil(length(Xe)/2));
    Xb1 = Xb(idx_b,:);
    Xb2 = Xb(setdiff(1:length(Xb),idx_b),:);
    Xe1 = Xe(idx_e,:);
    Xe2 = Xe(setdiff(1:length(Xe),idx_e),:);
    X1 = [Xb1;Xe1];
    X2 = [Xb2;Xe2];
    coef1 = get_coef(X1);
    if isempty(coef1)
        continue;
    end
    coef2 = get_coef(X2);
    if isempty(coef2)
        continue;
    end
    l = l+1;
    LRratio(l) = LRatio(X1,coef0,coef1)+LRatio(X2,coef0,coef2);
end
bd = prctile(LRratio,10);
if LRatio(Xb,coef0,coefb)+LRatio(Xe,coef0,coefe)<bd
    rej = true;
end
end

