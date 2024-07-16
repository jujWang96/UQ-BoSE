function []  = make_boxplot_bulleye(SNR,c,base,path)
if nargin==3
    path  = '~/Desktop/QE/simple_simulation/Bulleye/simulatedata_fix/';
end
m=4;
RIs = cell(1,3);
gRIs = cell(1,3);
source_mat = [];

for i = 1:length(c)
    rand_RIs = [];
    greedy_RIs = [];
    source_arr = [];
    for j = 1:length(base)
        source_arr = [source_arr;calc_num_in_source([1,3],0.03,c(i),base(j),SNR,2)];
        try
            if base(j)==100
                outputfile = strcat(path,'X_bulleye_simulated_ss1_grid_0.12_base',num2str(base(j)),'c',num2str(c(i)),'SNR',num2str(SNR),'m',num2str(m),'_seed');
            elseif base(j)==500
                outputfile = strcat(path,'X_bulleye_simulated_ss3_grid_0.1_base',num2str(base(j)),'c',num2str(c(i)),'SNR',num2str(SNR),'m',num2str(m),'_seed');

            else
                outputfile = strcat(path,'X_bulleye_simulated_ss5_grid_0.08_base',num2str(base(j)),'c',num2str(c(i)),'SNR',num2str(SNR),'m',num2str(m),'rand40_seed');
            end
        for l = 1:500
            load( strcat(outputfile,num2str(l),'.mat'))
            %greedy_reg_num(l) = length(find(~cellfun(@isempty,glb_selected_greedy)));
            greedy_RI(l) = RI_greedy;
            %random_reg_num(l) = length(find(~cellfun(@isempty,glb_selected)));
            rand_RI(l) = RI;
           
        
        end
        if isempty(rand_RIs)
            rand_RIs = rand_RI';  
            greedy_RIs = greedy_RI';
        else
            rand_RIs = [rand_RIs,rand_RI'];
            greedy_RIs = [greedy_RIs,greedy_RI' ];
        end
    %RIs = [RIs;greedy_RIs';rand_RIs'];
        catch
            disp(base(j))
            disp(c(i))
            disp(l)
            rand_RIs = [rand_RIs,nan(l,1)'];
            greedy_RIs = [greedy_RIs, nan(l,1)'];
        end
    end
    RIs{i} = rand_RIs;
    gRIs{i}= greedy_RIs;
    source_mat = [source_mat,source_arr];

end
l=500;
g = [repmat({'1'},l,1); repmat({'5'},l,1);repmat({'25'},l,1)];
sizeset = {strcat('\alpha=1, \beta','=', num2str(sum(source_mat(1,1:2)))), ...
    strcat('\alpha=5, \beta','=', num2str(sum(source_mat(2,1:2)))),...
    strcat('\alpha=25, \beta','=', num2str(sum(source_mat(3,1:2))))};
contrastset = {'c=1.4','c=2','c=2.8'};
boxplotGroup(RIs,'PrimaryLabels',contrastset,'SecondaryLabels',sizeset,'Colors',[0 .5 0; 0 0 1; .9 0 .9],'GroupType','betweenGroups','groupLines',true)
ylim([-0.2 1])