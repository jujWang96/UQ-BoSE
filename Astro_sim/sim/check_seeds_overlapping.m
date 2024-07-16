function violate = check_seeds_overlapping(seeds,seeds_pt,num_s,num_s_pt)
select_pt=[];
for i = 1:num_s
    select_pt = [select_pt,seeds{i}];
end
for i = 1:num_s_pt
    select_pt = [select_pt,seeds_pt{i}];
end 


if length(select_pt)>length(unique(select_pt))
    violate =1;
else
    violate=0;
end

