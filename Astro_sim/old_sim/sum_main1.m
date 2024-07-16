load('main1_rep_lambda_1000_num_200.mat')
tabulate(n_source)
mean(double(metric_all_frame), 1)
std(double(metric_all_frame), 0, 1)

load('main1_rep_lambda_1000_num_400.mat')
tabulate(n_source)
mean(double(metric_all_frame), 1)
std(double(metric_all_frame), 0, 1)

load('main1_rep_lambda_1000_num_600.mat')
tabulate(n_source)
mean(double(metric_all_frame), 1)
std(double(metric_all_frame), 0, 1)

load('main1_rep_lambda_2000_num_400.mat')
tabulate(n_source)
mean(double(metric_all_frame), 1)
std(double(metric_all_frame), 0, 1)

load('main1_rep_lambda_3000_num_600.mat')
tabulate(n_source)
mean(double(metric_all_frame), 1)
std(double(metric_all_frame), 0, 1)