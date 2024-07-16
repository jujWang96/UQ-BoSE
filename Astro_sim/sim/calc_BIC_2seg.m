function val = calc_BIC_2seg(X,target_curve,m)
%calculate the BIC value of the selected segmentation when there are only
%two segments 

log_like = -sum(log(1:length(X)))-length(X);
A = polyarea(target_curve(:,1),target_curve(:,2));
in = inpolygon(X(:,1),X(:,2),target_curve(:,1),target_curve(:,2));
out = ~in;
val = -2*(log_like+sum(in)*log(sum(in)/A)+sum(out)*log(sum(out)/(1-A)))+ m*log(length(X));



