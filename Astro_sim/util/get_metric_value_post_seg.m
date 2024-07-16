function val = get_metric_value_post_seg(n,num_region_nonempty, cell_area, cell_log_intensity, selected_nonempty, adj_mat)
[region_intensity, region_area, region_num_cells, adj_mat_region] =...
    get_region_int_connect(num_region_nonempty, cell_area, cell_log_intensity, selected_nonempty, adj_mat);
val = get_metric_value_flexible('log-like', n, selected_nonempty, region_area, region_intensity, region_num_cells);

