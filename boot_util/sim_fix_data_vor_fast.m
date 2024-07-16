function[new_X_set,Regions] = sim_fix_data_vor_fast(X,cell_area, n, DT,model,B,seed)
%generate bootstrapped date using voronoi tessellation
%input: 
%
% X: data point
% of the Voronoi vertices in V(has to be bounded using VoronoiLimit.m)
% model: 
%     |1. uniformly random
%     |2. Inhomogeneous possion
% B: number of bootstrap replicants. 
% seed: random seed.
%
%
%
%output: 
% new_X_set: a length-B cell array, with eahc cell containing a n by 2 matrix 
% of one bootstrap replicant. 
%
% The uniform points generating algorithm refers to ...
%https://math.stackexchange.com/questions/18686/uniform-random-point-in-triangle
%
rng(seed);

new_X_set = cell(B,1);
M = 1; %change M only for verification

[V, R] = voronoiDiagram(DT);

triagl = cell(n,1);
probs = cell(n,1);
Centers = zeros(n,2);
Regions = cell(B,1);
%partion voronoi cell into triangle and calculate each area
for i = 1:n
    v = V(R{i}, :);
    if isnan(cell_area(i))
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

infarea = find(isnan(cell_area));


%uniformly random generate
if model == 1 
    for b = 1:B
        N = n-length(infarea);%poissrnd(n);
        new_X = zeros(N*M,2);
        Regions = setdiff(1:n,infarea);
        for i = 1:N
            r = Regions(i);
            v = [V(R{r}, :);Centers(r,:)];
            p = probs{r};
            tri = triagl{r};
            %generate one random new point uniformly
            smplx = randsample(length(p),M,true,p);
            r1 = rand(M,1);
            pt = v(tri(smplx,1),:).*r1 + v(tri(smplx,2),:).*(1-r1);
            r2 = sqrt(rand(M,1));
            new_X((i-1)*M+1:i*M,:) = pt.*r2 + v(tri(smplx,3),:).*(1-r2); 
        end
        new_X = [new_X;X(infarea,:)];

        new_X_set{b} = double(new_X);
    end
    

% inhomogeneous poisson
elseif model == 2 
    for b = 1:B
        N = n-length(infarea);%poissrnd(n);
        new_X = zeros(N*M,2);
        Region = randsample(setdiff(1:n,infarea), N, true);
        Regions{b}= Region;
        for i = 1:N
            r = Region(i);                    
            v = [V(R{r},:);Centers(r,:)];
            p = probs{r};
            tri = triagl{r};
            %uniformly random sample in the selected area
            smplx = randsample(length(p),M,true,p);
            r1 = rand(M,1);
            pt = v(tri(smplx,1),:).*r1 + v(tri(smplx,2),:).*(1-r1);
            r2 = sqrt(rand(M,1));
            new_X((i-1)*M+1:i*M,:) = pt.*r2 + v(tri(smplx,3),:).*(1-r2); 
        end
        new_X = [new_X;X(infarea,:)];
        new_X_set{b} = double(new_X);
    end
    
% inhomogeneous poisson with original observations 
% if a cell is picked m times, m photons will be randomly placed in the
% cell; if if a cell is picked once, original photon will be placed 
elseif model == 3
    for b = 1:B
        N = n-length(infarea);%poissrnd(n);
        new_X = zeros(N*M,2);
        Region = randsample(setdiff(1:n,infarea), N, true);
        Regions{b}= Region;
        Region_id = unique(Region);
        for i = 1:length(Region_id)
           counts(i) = sum(Region==Region_id(i)); 
        end
        region_once = Region_id(counts==1);
        for i = 1:N
            r = Region(i);   
            if any(r == region_once)
                %the cell is selected only once, put the original X in 
                new_X((i-1)*M+1:i*M,:) = repmat([X(r,1),X(r,2)],M,1);
            else
                v = [V(R{r},:);Centers(r,:)];
                p = probs{r};
                tri = triagl{r};
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
    
            
else
    fprintf("Invalid model.");
end

    
        
end




