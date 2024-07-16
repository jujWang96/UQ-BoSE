function X = sim_inhomo_Pois_Gauss(range_x, range_y, lambda, loc, sigma, num, seed)
% simulate inhomogeneous Poisson process with several extended sources
% consisting of points/photons from a spatial Gaussian distribution
%
% Input variables:
%
% range_x: range of x
% range_y: range of y
% lambda: density of background
% loc: location of sources
% sigma: std of Gaussian distributions
% num: number of points/photons in each source
% seed: random seed
%
% Output variables:
%
% X: an n-by-2 matrix
%
% Examples:
%
% X = sim_inhomo_Pois_Gauss([0 1], [0 1], 100, [0.3 0.3; 0.7 0.7], [0.05 0.05], [100 100]);

if nargin==7
    rng(seed)
end

X = sim_homo_Pois(range_x, range_y, lambda);

% number of sources
n_s = length(sigma);

max_iter = 1e3;
for i = 1:n_s
    rand_num = mvnrnd(loc(i, :), sigma(i)^2*eye(2), max_iter);
    % make sure that the points are in the region
    index = find(rand_num(:, 1)>=range_x(1) & rand_num(:, 1)<=range_x(2) &...
        rand_num(:, 2)>=range_y(1) & rand_num(:, 2)<=range_y(2));
    index = index(1:num(i));
    rand_num = rand_num(index, :);
    % update X
    X = [X; rand_num];
end

figure
scatter(X(:, 1), X(:, 2), '.')
axis square

end