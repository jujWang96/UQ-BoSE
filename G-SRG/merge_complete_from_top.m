function [ log_like,region_intensity,region_area,selected] =...
    merge_complete_from_top(num, cell_area, ...
        cell_log_intensity, region_sets, adj_mat, n)
    %%
    %  merge subsets further to 2 segments so that
    %  likelihood is maximized among all merges
    %%

    [region_intensity, region_area, region_num_cells, adj_mat_region] =...
        get_region_int_connect(num, cell_area, cell_log_intensity, region_sets, adj_mat);
    %get all possible merge
    max_log = -realmax;
    max_seg = {};
    %search first half among all combinations
    for i = 1:floor(num/2)
        combos = combntns(1:num,i);
        for j = 1:length(combos)
            combon = setdiff(1:num,combos(j,:));
            if is_connect(combos(j,:)) && is_connect(combon)
                num1 = sum(region_intensity(combos(j,:)).*region_area(combos(j,:)));
                num2 = sum(region_intensity(combon).*region_area(combon));
                curr_log = num1*log(num1/sum(region_area(combos(j,:))))+...
                    num2*log(num2/sum(region_area(combon)));
                if curr_log>max_log
                    max_seg = {combos(j,:),combon};
                end
            end
        end

    end
    %update the regions after merge 
    if ~isempty(max_seg)
            %region_sets
       region_sets([max_seg{1}(2:end),max_seg{2}(2:end)])={0};      
         region_intensity(max_seg{1}(1)) = sum(region_intensity(max_seg{1}).*...
             region_area(max_seg{1}))/sum(region_area(max_seg{1}));
         region_intensity(max_seg{2}(1)) = sum(region_intensity(max_seg{2}).*...
             region_area(max_seg{2}))/sum(region_area(max_seg{2}));
         region_intensity([max_seg{1}(2:end),max_seg{2}(2:end)]) = 0;
         region_num_cells(max_seg{1}(1)) = sum(region_num_cells(max_seg{1}));
         region_num_cells(max_seg{2}(1)) = sum(region_num_cells(max_seg{2}));
         region_num_cells([max_seg{1}(2:end),max_seg{2}(2:end)]) = 0; 
         region_area(max_seg{1}(1)) = sum(region_area(max_seg{1}));
         region_area(max_seg{2}(1)) = sum(region_area(max_seg{2}));
         region_area([max_seg{1}(2:end),max_seg{2}(2:end)]) = 0;       
    end
     % compute the log-likelihood of the merged segmentation

    log_like = get_metric_value('log-like', n, region_sets, region_area, region_intensity, region_num_cells);
    selected = region_sets;
    %check if the subsets are connected 
    function res = is_connect(list)
        res = false;
        if length(unique(conncomp(graph(adj_mat_region(list,list)))))==1
            res = true;
        end

    end


    end