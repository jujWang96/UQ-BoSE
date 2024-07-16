function adj_mat = get_adj_mat_segment( seg_set, n )
% Gets adjacent matrix of a segmentation.
% 
% Args:
%   seg_set: A cell each element of which represents the a subgraph.
%   n: Total number of regions.
%
% Returns:
%   adj_mat: The adjacent matrix, where 1 represents that there is an edge,
%   and 0 represents that there is no edge.

adj_mat = false(n);
for i = 1:length(seg_set)
    for j = (i+1):length(seg_set)
        adj_mat(seg_set{i}, seg_set{j}) = true;
        adj_mat(seg_set{j}, seg_set{i}) = true;

    end
end

end