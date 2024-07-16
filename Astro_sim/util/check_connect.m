function flag = check_connect(k1, k2, region_sets, adj_mat)
% check connectivity between two region sets
%
% Input variables:
%
% k1, k2: indices of two region sets
% region_sets: a cell containing all region sets
% adj_mat: the adjacent matrix of points/voronoi cells
%
% Output variables:
%
% flag: whether the two region sets are connected

flag = 0;
for i = region_sets{k1}
    for j = region_sets{k2}
        % if there is one voronoi cell from each region set such that they
        % are connected
        if adj_mat(i, j)
            flag = 1;
            % exit the function
            return
        end
    end
end

end