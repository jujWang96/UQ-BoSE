function [bootstrap_X_set,N_in_out]= sim_fix_data(obs_X,range_x, range_y,polyset,B,seed)
%generate a dataset of fixed points with given segmentation and pdf
%Input: 
% X: original observations,N-by-2 matrix
% S: 1 by n cells of polygons
% p: probability of the point falling inside the segmentation.
% N: numebr of photons
% seed: random seed 
%Output:

if nargin==7
    rng(seed)
end
bootstrap_X_set = cell(B,1);
N_in_out = zeros(2,B);
[N,~] = size(obs_X);
%calculate the probability of following into each region 
for d = 1:length(polyset)
     num_photon(d) = sum(isinterior(pgon,obs_X(:, 1),obs_X(:, 2)));   
end
num_photon(d+1) = N-sum(num_photon);
 p = num_photon/N;
for b = 1:B
    innum = binornd( N , p );
    max_iter = N*10;

    % len_x = range_x(2)-range_x(1);
    len_x = range_x(2)-range_x(1);

    len_y = range_y(2)-range_y(1);
    %n represents the number of points, which follows a Poisson distribution
    % X is an n-by-2 matrix, where each row represents the location of one of
    % the points
    X(:, 1) = rand(max_iter, 1)*len_x+range_x(1);
    X(:, 2) = rand(max_iter, 1)*len_y+range_y(1);
    in = isinterior(X(:, 1),X(:, 2),polyset{1});
    Xin = X(in,:);
    if length(Xin)>=innum
        Xin = Xin(1:innum,:);
    else
        Xin = [Xin;sim_homo_fix(range_x, range_y, polyset{1}, innum-length(Xin), true)];
    end
    Xout = X(~in,:);
    if length(Xout)>=(N-innum)
        Xout = Xout(1:(N-innum),:);
    else
        Xout = [Xout; sim_homo_fix(range_x, range_y, polyset{2}, N-innum-lenth(Xout), true)];
    end
    N_in_out(:,b) = [length(Xin);length(Xout)];
   
    bootstrap_X_set{b} = [Xin; Xout];
end
end
%---------------------------------------------------------
function X = sim_homo_fix(range_x, range_y, pgon, n, in)
% simulate homogeneous uniform distribution with fixed n in given range 
% Input variables:
%
% range_x: the range of x
% range_y: the range of y
% S:L-by-2 segmentation contour.
% in:boolean variable, if true then generate data inside S, otherwise
% outside 
%
% Output variables:
%
% X: an n-by-2 matrix, where each row represents the location of one of the 
% points

[x,y] = boundary(pgon);
if in 
    lenx = max(x)-min(x);
    leny = max(y)-min(y);
    X(:, 1) = rand(n*10, 1)*lenx+min(x);
    X(:, 2) = rand(n*10, 1)*leny+min(y);
    X = X(isinterior(pgon,X(:, 1),X(:, 2)),:);
    X = X(1:n,:);
else
    lenx = range_x(2)-range_x(1);

    leny = range_y(2)-range_y(1);
    X(:, 1) = rand(n*10, 1)*lenx+range_x(1);
    X(:, 2) = rand(n*10, 1)*leny+range_y(1);
    X = X(~isinterior(pgon,X(:, 1),X(:, 2)),:);
    X = X(1:n,:);
end


end