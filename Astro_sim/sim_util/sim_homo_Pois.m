%author: Minjie Fan@UCDavis/Google

function X = sim_homo_Pois(range_x, range_y, lambda, seed)
% simulate homogeneous Poisson process
% Input variables:
%
% range_x: the range of x
% range_y: the range of y
% lambda: density of the process (which is a constant)
% seed: random seed
%
% Output variables:
%
% X: an n-by-2 matrix, where each row represents the location of one of the 
% points

len_x = range_x(2)-range_x(1);
len_y = range_y(2)-range_y(1);

if nargin==4
    rng(seed)
end
% n represents the number of points, which follows a Poisson distribution
n = poissrnd(len_x*len_y*lambda);
% X is an n-by-2 matrix, where each row represents the location of one of
% the points
X = zeros(n, 2);
X(:, 1) = rand(n, 1)*len_x+range_x(1);
X(:, 2) = rand(n, 1)*len_y+range_y(1);

% scatter(X(:, 1), X(:,2), '.')
% axis square

end