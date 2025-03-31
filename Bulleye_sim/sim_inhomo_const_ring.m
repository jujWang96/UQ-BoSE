
function X = sim_inhomo_const_ring(range_x, range_y, loc, radius, alpha,lambda,num_in_source, seed, show_plot)
% simulate an inhomogeneous Poisson process with several homogeneous, 
% round-shaped extended sources with rings around
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
% X = sim_inhomo_const_ring([0 1], [0 1], 100, [0.3 0.3; 0.7 0.7], [0.1 0.05; 0.1,0.02], [100 10; 100 30]);

if nargin==7
    show_plot = false;
elseif nargin==8
    rng(seed)
    show_plot = false;
elseif nargin==9
    rng(seed)
end

[n_s,D] = size(radius);
totalN = (sum(num_in_source)+lambda)*alpha;

pois_num = mnrnd(totalN,double([num_in_source,lambda]/sum([num_in_source,lambda])));


X = sim_homo_const(range_x, range_y, pois_num(end));
 
% number of sources
num = pois_num(1:(end-1));
for i = 1:n_s
    max_iter = max(num(i,:))*10;
    rand_circle = rand(max_iter, 2);
    % generate random points within the square
    rand_circle(:, 1) = rand_circle(:, 1)*radius(i,1)*2+loc(i,1)-radius(i,1);
    rand_circle(:, 2) = rand_circle(:, 2)*radius(i,1)*2+loc(i,2)-radius(i,1);
    % find the index of the points that are located in the circle
    index = find(((rand_circle(:, 1)-loc(i, 1)).^2+(rand_circle(:, 2)-loc(i, 2)).^2)<=radius(i)^2);
    index = index(1:num(i,1));
    rand_circle = rand_circle(index, :);
    % update X
    X = [X; rand_circle];
    in_rad = radius(i,1);
    for d =2:D
        if radius(i,d)==0
            break
        end
        rand_ring = unif_rand_ring(num(i,d),loc(i,:),in_rad,radius(i,d));

        X = [X; rand_ring];
        in_rad = in_rad+radius(i,d);
    end
    
end

if show_plot
    figure
    scatter(X(:, 1), X(:, 2), 1)
    axis square
end
    function rand_ring = unif_rand_ring(num_ring, loc,in_rad,width)
        rho = sqrt(rand(num_ring,1)*((width+in_rad)*(width+in_rad)-in_rad*in_rad)+in_rad*in_rad) ;
        theta = rand(num_ring,1)*2*pi;
        rand_ring = [rho.*cos(theta)+loc(1),rho.*sin(theta)+loc(2)];
        
    end
end
