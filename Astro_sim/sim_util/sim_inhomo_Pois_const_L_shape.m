function X = sim_inhomo_Pois_const_L_shape(range_x, range_y, num, seed, show_plot)
% simulate an inhomogeneous Poisson process with several homogeneous, 
% round-shaped extended sources
% 
% Input variables:
%
% range_x: range of x
% range_y: range of y
% lambda: density of background
% loc: location of sources
% radius: radius of sources
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

if nargin==3
    show_plot = false;
elseif nargin==4
    rng(seed)
    show_plot = false;
elseif nargin==5
    rng(seed)
end

% number of sources
n_s = size(range_x, 1);

X = [];
len_x = range_x(:, 2)-range_x(:, 1);
len_y = range_y(:, 2)-range_y(:, 1);

for i = 1:n_s
    rand_num = zeros(num(i), 2);
    % generate random points within the square
    rand_num(:, 1) = rand(num(i), 1)*len_x(i)+range_x(i, 1);
    rand_num(:, 2) = rand(num(i), 1)*len_y(i)+range_y(i, 1);
    % update X
    X = [X; rand_num];
end

if show_plot
    figure
    scatter(X(:, 1), X(:, 2), '.')
    axis square
end

end