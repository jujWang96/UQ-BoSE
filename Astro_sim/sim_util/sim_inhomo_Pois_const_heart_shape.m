function X = sim_inhomo_Pois_const_heart_shape(range_x, range_y, lambda, loc, radius, num, seed, show_plot)
% simulate an inhomogeneous Poisson process with several homogeneous, 
% heart-shaped extended sources
% 
% Input variables:
%
% range_x: range of x
% range_y: range of y
% lambda: density of background
% loc: location of sources
% radius: radius of sources, radius have to be less than 0.2 for the shape
% to be compatible with the region
% num: number of points/photons in each source
% seed: random seed
% show_plot: whether show plot
%
% Output variables:
%
% X: an n-by-2 matrix
%
% Examples:
%
% X = sim_inhomo_Pois_const([0 1], [0 1], 100, [0.3 0.3; 0.7 0.7], [0.1 0.1], [100 100]);

if nargin==6
    show_plot = false;
elseif nargin==7
    rng(seed)
    show_plot = false;
elseif nargin==8
    rng(seed)
end

X = sim_homo_Pois(range_x, range_y, lambda);

% number of sources
n_s = length(radius);

max_iter = max(num)*10;
for i = 1:n_s
    inum = poissrnd(num(i));
    rand_num = rand(max_iter, 2)-0.5;
    % generate random points within the square
    [rand_theta,rand_rho] = cart2pol(rand_num(:, 1),rand_num(:, 2));
    index = find(rand_rho<=2*radius*(1-cos(rand_theta)));

    rand_num(:, 1) = rand_num(:, 1)+loc(i, 1);
    rand_num(:, 2) = rand_num(:, 2)+loc(i, 2);
    % find the index of the points that are located in the circle
    %index = find(((rand_num(:, 1)-loc(i, 1)).^2+(rand_num(:, 2)-loc(i, 2)).^2-1).^3<=((rand_num(:, 1)-loc(i, 1)).^2.*(rand_num(:, 2)-loc(i, 2)).^3));
    index = index(1:inum);
    rand_num = rand_num(index, :);
    % update X
    X = [ X;rand_num];
end

if show_plot
    figure
    scatter(X(:, 1), X(:, 2), '.')
    axis square
    axis( [0 1 0 1]);
end

end