function X = sim_inhomo_const_general(range_x, range_y, lambda, x,y, num, seed, show_plot)
% simulate an inhomogeneous Poisson process with the given object
% boundaries 
% Input variables:
%
% range_x: range of x
% range_y: range of y
% lambda: density of background
% x: cells of boundry coordinates of x
% y: cells of boundry coordinates of y

% num: number of points/photons in each source
% seed: random seed
% show_plot: whether show plot
%
% Output variables:
%
% X: an n-by-2 matrix
%


if nargin==6
    show_plot = false;
    X = sim_homo_const(range_x, range_y, lambda);

elseif nargin==7
    rng(seed)
    show_plot = false;
    X = sim_homo_const(range_x, range_y, lambda,seed);

elseif nargin==8
    rng(seed)
    
    X = sim_homo_const(range_x, range_y, lambda,seed);

end
length(X);
% number of sources
n_s = length(x);

max_iter = max(num)*10;
for i = 1:n_s
    pgon = polyshape(x{i},y{i});
    rand_num = rand(max_iter, 2);
    % generate random points within the rectangle of the polygon 
    rand_num(:, 1) = rand_num(:, 1)*(max(x{i})-min(x{i}))+min(x{i});
    rand_num(:, 2) = rand_num(:, 2)*(max(y{i})-min(y{i}))+min(y{i});
    
    % find the index of the points that are located in the circle
    index = isinterior(pgon,rand_num(:, 1),rand_num(:, 2));
    rand_num = rand_num(index, :);
    rand_num = rand_num(1:num(i),:);
    % update X
    X = [X; rand_num];
end

if show_plot
    figure
    scatter(X(:, 1), X(:, 2), '.')
    axis square
end

end