function [init_sets, labeled_set] = SRG_graph( init_sets, flux, area_cell, n, adj_mat, invalid)

% graph-based seeded region growing
% init_set: the initial growing region sets
% flux: the fluxes of all voronoi cells
% n: the number of voronoi cells

% initialize all the growing regions
m = length(init_sets);
region_flux = zeros(1, m);
region_area = zeros(1, m);
labeled_set = zeros(1, m);
nbr_set = cell(m, 1);
for j = 1:m
    index = init_sets{j}(1);
    labeled_set(j) = index;
    region_flux(j) = flux(index);
    region_area(j) = area_cell(index);
end

for j = 1:m
    nbr_set{j} = cmp_nbr( init_sets{j}, adj_mat, n, [labeled_set invalid] );
end

for i = 1:n-m-length(invalid)
    
    if mod(i, 100)==0
        i
    end
    
    % select the pair of a growing region and a cell with smallest value of
    % the criterion
    min_c = realmax;
    for j = 1:m
        for k = nbr_set{j}
            c = abs(flux(k)-region_flux(j));
            if c<min_c
                min_c = c;
                min_j = j;
                min_k = k;
            end
        end
    end
    
    % update the selected region
    labeled_set = [labeled_set min_k];
    init_sets{min_j} = [init_sets{min_j} min_k];
    region_flux(min_j) = log((exp(region_flux(min_j))*region_area(min_j)+exp(flux(min_k))*...
        area_cell(min_k))/(region_area(min_j)+area_cell(min_k)));
    region_area(min_j) = region_area(min_j)+area_cell(min_k);
    % add the neighbors of min_k to min_j
    nbr_set{min_j} = [nbr_set{min_j} cmp_nbr( min_k, adj_mat, n, [labeled_set invalid] )];
    % delete min_k from neighbor sets
    for j = 1:m
        nbr_set{j} = setdiff(nbr_set{j}, min_k);
    end
    
end
        
end