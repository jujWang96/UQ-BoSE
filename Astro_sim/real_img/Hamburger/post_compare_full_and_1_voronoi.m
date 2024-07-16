clear 
close all
filename = 'real_full_to_1_random_voronoi';
fid=fopen('Hamburger_11007/hamburger_0123.txt');
C=textscan(fid,'%f%f%f%f%f%f%f%f%f','HeaderLines',11,'delimiter','\r\n','CommentStyle','*');
fclose(fid);
penalty_term = 6;
M = 150; %run greedy merge until 150 region left
X_full = [C{1},C{2}];
frame = C{8};
arcsec_full = [C{4},C{5}];
fid=fopen('Hamburger_11007/hamburger_1.txt');
C=textscan(fid,'%f%f%f%f%f%f%f%f%f','HeaderLines',11,'delimiter','\r\n','CommentStyle','*');
fclose(fid);

X1 = [C{1},C{2}];
arcsec = [C{4},C{5}];
%invalid_full = ismember(arcsec_full(:,1),arcsec(invalid,1)) & ismember(arcsec_full(:,2),arcsec(invalid,2));
load('real_1_voronoi_result_penalty6_2023_2_16_20_25_29.mat')
invalid_1 = invalid;
cell_log_intensity_1 = cell_log_intensity;
cell_area_1 = cell_area;
%cell_log_intensity_1 = cell_log_intensity(valid);
load('real_full_voronoi_result_penalty6_2023_2_12_22_36_53.mat')
selected_full_convert = cell([]);
idx = 0;
cell_log_intensity_1_valid = [];
cell_log_intensity_full_convert_valid = [];
cell_area_1_valid = [];
fig = figure;
hold on
for i = 1:num
    if ~isempty(selected{i})
        seg_in_frame = ismember(arcsec(:,1), arcsec_full(selected{i},1)) & ismember(arcsec(:,2), arcsec_full(selected{i},2));
        seg_ids = find(seg_in_frame==1);
        seg_ids = setdiff(seg_ids,invalid_1);
        if ~isempty(seg_ids)
            idx = idx+1;

            selected_full_convert{idx} = seg_ids;
            seg_in_frame_full = ismember(arcsec_full(:,1),arcsec(seg_ids,1)) & ismember( arcsec_full(:,2),arcsec(seg_ids,2));
            seg_ids_revert = find(seg_in_frame_full==1);
            %scatter(X_full(seg_ids_revert,1),X_full(seg_ids_revert,2))
            cell_area_1_valid = [cell_area_1_valid; cell_area_1(seg_ids)];
            cell_log_intensity_1_valid = [cell_log_intensity_1_valid; cell_log_intensity_1(seg_ids)];
            biased = find(cell_log_intensity_1(seg_ids)-cell_log_intensity(seg_ids_revert)<-1.5);
            scatter(X1(seg_ids(biased),1),X1(seg_ids(biased),2),'k')
            cell_log_intensity_full_convert_valid = [cell_log_intensity_full_convert_valid; cell_log_intensity(seg_ids_revert)];
        end
    end
end
axis equal
saveas(fig, strcat(filename,'_log_int_different'), 'epsc')

load('real_1_voronoi_result_penalty6_2023_2_16_20_25_29.mat')
figure
plot_segmentation_wo_voronoi(DT, selected_full_convert, cx, cy,colors ,8,true)
fig = figure;
triplot(DT, 'Color', GRAY)
hold on
for i = 1:length(selected_full_convert)
    log_int = log(sum(exp(cell_log_intensity(selected_full_convert{i})).*cell_area(selected_full_convert{i}))/sum(cell_area(selected_full_convert{i})));
    scatter(cx(selected_full_convert{i}), cy(selected_full_convert{i}), 12, log_int*ones(length(selected_full_convert{i}), 1), 'filled')
end
colorbar('EastOutside')
colormap(hsv)
caxis([9 16])
axis image
set(gca, 'fontsize', 12)
min_white_margin(gca);

saveas(fig, strcat(strcat(filename,'_penalty',num2str(penalty_term)),'segmentat_result_region'), 'epsc')

num_region_convert = length(selected_full_convert);
val = get_metric_value_post_seg(n,num_region_convert, cell_area, cell_log_intensity, selected_full_convert, adj_mat);
BIC_convert = -2*val+penalty_term*(num_region_convert-1)*log(n);
index = 0;
selected_nonempty = {};
for i = 1:length(selected)
    if ~isempty(selected{i})
        index = index + 1;
        selected_nonempty{index} = selected{i};
    end
end

num_region_origin = length(selected_nonempty);
val = get_metric_value_post_seg(n,num_region_origin, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
BIC_origin = -2*val+penalty_term*(num_region_origin-1)*log(n);
%plot the voronoi cells intensity 
figure
ecdf(cell_log_intensity_1+log(4))
hold on
ecdf(cell_log_intensity_full_convert_valid)

fig = figure
histogram(cell_log_intensity_1_valid-cell_log_intensity_full_convert_valid,50);
saveas(fig, strcat(filename,'hist_log_int'),'epsc')