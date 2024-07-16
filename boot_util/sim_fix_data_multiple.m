function [bootstrap_X_set]= sim_fix_data_multiple(obs_X,range_x, range_y,polyset,B,seed,N)
%generate a dataset of fixed number of point with given segmentation and pdf
%Input: 
% obs_X: original observations
% polyset: the polygon representing the segmentation of ROI
% N: total numebr of photons
% seed: random seed 
%Output:
% X:N-by-2 matrix


% [N,~] = size(obs_X)
if nargin==7
    rng(seed)
elseif nargin==6
   rng(seed)
   [N,~] = size(obs_X);
else
   [N,~] = size(obs_X);
end
bootstrap_X_set = cell(B,1);
for d = 1:length(polyset)
     num_photon_estimate(d) = sum(isinterior(polyset{d},obs_X(:, 1),obs_X(:, 2)));   
end
p = num_photon_estimate/N;
for b = 1:B
    X_in = [];
    num_photon = mnrnd(N,p,1);
    %generate photons in each segments
    for d = 1:length(num_photon)
       
        X_in =[X_in; sim_homo_fix(range_x, range_y, polyset{d}, num_photon(d), true)];
    end
    bootstrap_X_set{b} = X_in;
    %disp(length(X_in))
end
end
%---------------------------------------------------------
function X = sim_homo_fix(range_x, range_y, pgon, n, in)
% simulate homogeneous uniform distribution with fixed n in given range 
% Input variables:
%
% range_x: the range of x
% range_y: the range of y
% pgon: polygon shape
% in:boolean variable, if true then generate data inside S, otherwise
% outside 
%
% Output variables:
%
% X: an n-by-2 matrix, where each row represents the location of one of the 
% points
X = [];
n_left = n; 
[x,y] = boundary(pgon);

if in 
    %generate photon inside polygon
    lenx = max(x(:))-min(x(:));
    leny = max(y(:))-min(y(:));
    while length(X)<n
        clear tempt_X
        tempt_X(:, 1) = rand(n_left*10, 1)*lenx+min(x(:));
        tempt_X(:, 2) = rand(n_left*10, 1)*leny+min(y(:));
        in = isinterior(pgon,tempt_X(:, 1),tempt_X(:, 2));
        tempt_X = tempt_X(in,:);
        if length(tempt_X)<n_left
            X = [X;tempt_X];
            n_left = n-length(X);
        else
            X = [X;tempt_X(1:n_left,:)];
        end
    end
%generate photon outside polygon
else
    while length(X)<n
        clear tempt_X

        lenx = range_x(2)-range_x(1);

        leny = range_y(2)-range_y(1);
        tempt_X(:, 1) = rand(n_left*10, 1)*lenx+range_x(1);
        tempt_X(:, 2) = rand(n_left*10, 1)*leny+range_y(1);
        in = isinterior(pgon,tempt_X(:, 1),tempt_X(:, 2));
        tempt_X = tempt_X(~in,:);
        if length(tempt_X)<n_left
            X = [X;tempt_X];
            n_left = n-length(X);
        else
            X = [X;tempt_X(1:n_left,:)];
        end
    end
end


end