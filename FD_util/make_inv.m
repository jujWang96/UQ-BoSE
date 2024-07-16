%function [coef_a,coef_b] = make_inv(coef,Mp)
function coef = make_inv(coef)

N = length(coef);
if nargin==1
    Mp = floor((N-1)/2);
end
coef = make_trans_inv(coef);
coef = make_scale_inv(coef);
coef = make_startpt_inv2(coef);
%[coef_a,coef_b] = make_startpt_inv(coef,Mp);
%coef_a = make_rotate_inv(coef_a);
%coef_b = make_rotate_inv(coef_b);
end
