function [new_X_set,Regions] = sim_fix_data_vor_fast_vk(X,cell_area, n, DT,B,seed)
%generate bootstrapped date using voronoi tessellation.
% For each bootstrap data, first sample n points with replacement, use the
% uniquely selected point to create new voronoi tessellation. If the point
% is selected once, put a observation in the original place; if a
% point is selected K times, put K observations in the corresponding
% voronoi cell.
%input: 
%
% X: data point
% of the Voronoi vertices in V(has to be bounded using VoronoiLimit.m)
% B: number of bootstrap replicants. 
% seed: random seed.
%
%output: 
% new_X_set: a length-B cell array, with each cell containing a n by 2 matrix 
% of one bootstrap replicant. 
 
% The uniform points generating algorithm refers to ...
%https://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle
%
rng(seed);

new_X_set = cell(B,1);
M = 1; %change M only for verification

[V, R] = voronoiDiagram(DT);
infarea = find(isnan(cell_area));

for b = 1:B
    N = n-length(infarea);%poissrnd(n);
    new_X = zeros(N*M,2);
    Region = randsample(setdiff(1:n,infarea), N, true);
    Regions{b}= Region;
    Region_id = unique(Region);
    counts = zeros(1,length(Region_id));
    for i = 1:length(Region_id)
       counts(i) = sum(Region==Region_id(i)); 
    end
    uniqueN = length(Region_id)+length(infarea);
    triagl = cell(uniqueN,1);
    probs = cell(uniqueN,1);
    Centers = zeros(uniqueN,2);
    [cxb, cyb, nb, DTb, Eb, cell_log_intensityb, cell_areab] = init_comp(double([X(Region_id,:);X(infarea,:)]), [0 1], [0 1], ones(size([X(Region_id,:);X(infarea,:)], 1), 1));

    [Vb, Rb] = voronoiDiagram(DTb);
       
    %partion voronoi cell into triangle and calculate each area
    for i = 1:uniqueN
        v = Vb(Rb{i}, :);
        if isnan(cell_areab(i))
            continue;
        end
        CH = convhulln(v);
        length_v = length(v);
        length_CH = length(CH);
        
        cent = mean(v,1);
        v(length_v+1,:) = cent;
        %save the partition to triagl
        tri = [CH,repmat(length_v+1,length_CH,1)];
        triagl{i} = tri;
        p = zeros(1,length_CH);
        for j = 1:length_CH
            p(j) = abs(det(v(tri(j,1:2),:) - cent));
        end
        p = p/sum(p);
        probs{i} = p;
        %save the center to V and update R 
        Centers(i,:) = cent;
    end
    region_once = Region_id(counts==1);
    
    for i = 1:N
        r = Region(i);   
        rnew =find(Region_id==r);
        if any(r == region_once) || isnan(cell_areab(rnew))
            %the cell is selected only once, put the original X in 
            new_X((i-1)*M+1:i*M,:) = repmat([X(r,1),X(r,2)],M,1);
        else
            v = [Vb(Rb{rnew},:);Centers(rnew,:)];
            p = probs{rnew};
            tri = triagl{rnew};
            %uniformly random sample in the selected area
            smplx = randsample(length(p),M,true,p);
            r1 = rand(M,1);
            pt = v(tri(smplx,1),:).*r1 + v(tri(smplx,2),:).*(1-r1);
            r2 = sqrt(rand(M,1));
            new_X((i-1)*M+1:i*M,:) = pt.*r2 + v(tri(smplx,3),:).*(1-r2);
        end
    end
    new_X = [new_X;X(infarea,:)];
    
    new_X_set{b} = double(new_X);
end

    
        
end




