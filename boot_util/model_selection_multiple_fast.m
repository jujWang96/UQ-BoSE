

function [smooth_polyset,min_pair_num] = model_selection_multiple_fast(X,target_x,target_y,structure,out_bd_x,out_bd_y,IC,M_pairs)
%perform the OPTIMIZED model selection with parallel for FDs of multiple target curves assuming that the target
%contours will not intersect with the inner bound or outer bound. 
%**could deal with model than one disconnected boundries in each stage
%input:
%   X: 2d photons
%   inner_bd, out_bd: 2d arrays of curve for inner and outer bound
%   target_x: 1-by-d cells of x coordinates
%   target_y: 1-by-d cells of y coordinates
%   structure: 1x(d-1) cell of the indices of x/y coordinates cells that
%   construct one connected polygon
%   IC: criterion for model selection
%smooth_curves=0;
warning('off')
[n,~] = size(X);
for d = 1:length(target_x)

    [centx,centy] = centroid(polyshape(target_x{d},target_y{d}));
    [x_sample,y_sample] = sample_curve(target_x{d},target_y{d},300,centx, false);
    f_coefs(:,d) = FD(x_sample,y_sample);
    
end

if nargin == 7 && strcmp(IC, 'BIC')
    M_pairs = ones(1,d)*2;
elseif nargin == 7 && strcmp(IC, 'AIC')
    M_pairs = ones(1,d)*5;

end
L = length(structure);
%two sets of length L that score the area and number of events in each
%possible
%smoothed layer 
area_set = cell(1,L);
num_set = cell(1,L);
tic
for l = 1:L
    structure_l = structure{l};

    layer_area = zeros(1,prod(M_pairs(structure_l)));
    layer_num = zeros(1,prod(M_pairs(structure_l)));
    M_pair = M_pairs(structure_l);

    parfor idx = 1:prod(M_pairs(structure_l))
        
        %constructure the smoothed polygon of cooresponding number of FD
        %pairs
        mat_idx = get_mat_index(idx,M_pair);
        smooth_x = cell(1,length(mat_idx));
        smooth_y = cell(1,length(mat_idx));
        for n_curve = 1:length(mat_idx)
            [smooth_x{n_curve},smooth_y{n_curve}]=iFD(f_coefs(:,structure_l(n_curve)),mat_idx(n_curve)*2+1);
        end
        if l<L     
            pgon = polyshape(smooth_x,smooth_y);
        else
            %this is the outer most layer
            pgon = polyshape([smooth_x,{out_bd_x}],[smooth_y,{out_bd_y}]);
        end
        
        %check if any curve intersect with the other 

        if pgon.NumRegions>1
            continue 
        end
        %store the corresponding area and number of events 
        layer_area(idx) = area(pgon);
        layer_num(idx) = sum(isinterior(pgon,X(:,1),X(:,2)));
   
    end
    area_set{l} = layer_area;
    num_set{l} = layer_num;
    disp(strcat('finish layer ',num2str(l)))
end

toc


tic
%calculate the loglikelihood value of all FD pairs combinations and keep track of the FD pairs returns
%smallest value for all required 
min_val = Inf*ones(1,length(IC));
min_pair_num = cell(1, length(IC));
for idx = 1:prod(M_pairs)
        log_like = -sum(log(1:n))-n;
        
        %decompose to the pair of each boundary
        mat_idx = get_mat_index(idx,M_pairs);
        for l =1:L
            arr_idx = get_arr_index(mat_idx(structure{l}),M_pairs(structure{l}));
            %skip the current combination of pairs
            if num_set{l}(arr_idx)==0
                log_like = -Inf;
                break
            end
            log_like = log_like+num_set{l}(arr_idx)*log(num_set{l}(arr_idx)/area_set{l}(arr_idx));
        end
        
        
        %update the value for all specific model selection criterion
        for IC_idx = 1:length(IC)
            if strcmp(IC(IC_idx), 'BIC')
                val = (-2)*log_like+log(n)*(2*(2*sum(mat_idx)+length(mat_idx))+1);
            elseif strcmp(IC(IC_idx), 'AIC')
                val = (-2)*log_like+2*(2*(2*sum(mat_idx)+length(mat_idx))+1);
            elseif strcmp(IC(IC_idx), 'MDL')
                val = (-2)*log_like;
                for l =1:L
                    arr_idx = get_arr_index(mat_idx(structure{l}),M_pairs(structure{l}));
                    %skip the current combination of pairs if not valid 
                    if num_set{l}(arr_idx)==0
                        val = Inf;
                        break
                    end
                    val = val+log(num_set{l}(arr_idx))*2*(2*sum(mat_idx(structure{l}))+length(structure{l}));
                end
            else
                disp('criterion not listed')
                val=Inf;
            end

            if val<min_val(IC_idx)
                min_val(IC_idx) = val;
                min_pair_num{IC_idx} = mat_idx;
            end
        end
        
        
end
toc
smooth_polyset = cell(1,length(IC));

%construct the smoothed polyset for each criterion

for IC_idx = 1:length(IC)
    if isempty(min_pair_num{IC_idx})
        continue
    end
    smooth_poly = cell(1,L);
    for l = 1:L
        clear smooth_x smooth_y 
        for n_curve = 1:length(structure{l})    
            [smooth_x{n_curve},smooth_y{n_curve}]=iFD(f_coefs(:,structure{l}(n_curve)),min_pair_num{IC_idx}(structure{l}(n_curve))*2+1);
            if l<L     
                smooth_poly{l} = polyshape(smooth_x,smooth_y);
            else
                %this is the outer most layer
                smooth_poly{l} = polyshape([smooth_x,{out_bd_x}],[smooth_y,{out_bd_y}]);
            end
        end

    end
    smooth_polyset{IC_idx} = smooth_poly;

end

   
end
 function mat_idx = get_mat_index(arr_idx,M_pairs)
    %given the 1-d array index and the size of matrix, return the index in
    %the multiple dimension array
        mat_idx = ones(1,numel(M_pairs));
        curr_dim = numel(M_pairs);
        while curr_dim>1
            if arr_idx>prod(M_pairs(1:(curr_dim-1)))
                if mod(arr_idx,prod(M_pairs(1:(curr_dim-1))))~=0
                    mat_idx(curr_dim) = mat_idx(curr_dim)+floor(arr_idx/prod(M_pairs(1:(curr_dim-1))));
                   
                    arr_idx = mod(arr_idx,prod(M_pairs(1:(curr_dim-1))));
                else
                    mat_idx(curr_dim) = mat_idx(curr_dim)+floor(arr_idx/prod(M_pairs(1:(curr_dim-1))))-1;
                    arr_idx=prod(M_pairs(1:(curr_dim-1)));
                end
            end
            curr_dim = curr_dim-1;  
        end
        mat_idx(curr_dim) = arr_idx;
    end
    function arr_idx = get_arr_index(mat_idx,M_pairs)
    %given the multiple dimension array index and the size of matrix, 
    %return the index in the 1-d array index
         arr_idx = 0;
         curr_dim = numel(M_pairs);
         while curr_dim>1
             arr_idx = arr_idx+(mat_idx(curr_dim)-1)*prod(M_pairs(1:(curr_dim-1))); 
             curr_dim = curr_dim-1;
         end
         arr_idx= arr_idx+mat_idx(1);  
    end
     
