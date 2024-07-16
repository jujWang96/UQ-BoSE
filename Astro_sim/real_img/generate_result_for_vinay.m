% load fitted result
load('real_full_result_2018_11_1_23_4_37.mat')

% Comments from Vinay: While playing with segment_id.txt, I realized that 
% it would be more useful with a couple of extra columns ? one is the 
% voronoi cell area associated with each photon (so we can get a quick 
% estimate of the average intensity in each segment), and the time, and 
% either the energy or PI value that comes with each photon. 
% read data (in its original scale)
% CIAO command:
% dmlist "./315/repro/acisf00315_repro_evt2.fits[cols x, y, TIME, PI][sky=region(ds9.reg)]" \
% data,clean outfile=photons
filename='photon_loc_time_pi.txt';
delimiterIn = ' ';
headerlinesIn = 1;
A = importdata(filename,delimiterIn,headerlinesIn);

X = A.data;
n_cols = size(X, 2);
clear A

% normalization factor
min_x = min(X(:, 1));
max_x = max(X(:, 1));
norm_factor = max_x - min_x;
cell_area_denorm = cell_area * norm_factor^2;
X(:, n_cols + 1) = cell_area_denorm; 

% init the last column (recording the segment id)
X(:, n_cols + 2) = 0;

num = length(selected);
area_all = [];
log_int_all = [];
selected_nonempty = {};
index = 0;

for i = 1:num
    if ~isempty(selected{i})
        index = index+1;
        selected_nonempty{index} = selected{i};
        area = sum(cell_area(selected{i}));
        area_all = [area_all, area];
        log_int = log(sum(exp(cell_log_intensity(selected{i})).*cell_area(selected{i}))/area);
        log_int_all = [log_int_all, log_int];
    end
end

% sorted by log intensity descendingly
[~, sorted_index] = sort(log_int_all, 'descend');
selected_nonempty = selected_nonempty(sorted_index);

% segment_id starts from 1. segment_id = 0 means the photon does not belong
% to any segments. They are near the edge.
num_segments = length(selected_nonempty);
for segment_id = 1:num_segments
    X(selected_nonempty{segment_id}(1, :), n_cols + 2) = segment_id;
end

X_table = array2table(X, 'VariableNames', ...
    {'x' 'y' 'time', 'PI', 'cell_area', 'segment_id'});
writetable(X_table, 'segment_info.txt', 'Delimiter', ' ')  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for segment_id = 1:num_segments
    h = figure;
    scatter(X(:, 1), X(:, 2), '.')
    hold on
    x_segment = X(X(:, n_cols + 2) == segment_id, 1);
    y_segment = X(X(:, n_cols + 2) == segment_id, 2);
    convhull_points = convhull(x_segment, y_segment);
    scatter(x_segment, y_segment, '.')
    plot(x_segment(convhull_points), y_segment(convhull_points), '-')
    axis tight
    title(['Segment ', num2str(segment_id)])
    saveas(h, ['convhull_segment_', num2str(segment_id), '.png'])
    close all
end

% write to movie
writerObj = VideoWriter('segments.avi');
writerObj.FrameRate = 1;
open(writerObj);
for segment_id = 1:num_segments
  thisimage = imread(['convhull_segment_', num2str(segment_id), '.png']);
  writeVideo(writerObj, thisimage);
end
close(writerObj);
