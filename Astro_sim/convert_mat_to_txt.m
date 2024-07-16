clear;
parpath = '~/src/gsrg/';
plotpath = 'Astro_sim/real_data/plots/';
resultpath = 'Astro_sim/real_data/results/';
addpath(genpath(strcat(parpath,'Astro_sim')))
addpath(genpath(strcat(parpath,'G-SRG')))
addpath(genpath(strcat(parpath,'plot_util')))
addpath(genpath(strcat(parpath,'contour_util')))
addpath(genpath(strcat(parpath,'FD_util')))
addpath(genpath(strcat(parpath,'boot_util')))
addpath(genpath(strcat(parpath,'gcr_util')))
addpath(genpath(strcat(parpath,'Arp299XMM')))
addpath(genpath(strcat(parpath,'NGC2300')))
addpath(genpath(strcat(parpath,'NGC2300XMM')))



load('Arp299_MOS1_evt_0.5-8.0keV_scaled.mat');
max(X)
min(X)
dlmwrite('Arp299_MOS1_evt_0.5-8.0keV.txt', X, 'delimiter', '\t');
load('ngc2300_box_058kev_evt2.mat');
max(X)
min(X)
dlmwrite('ngc2300_box_058kev_evt2.txt', X, 'delimiter', '\t');
load('ngc2300_MOS1_evt_0.5-8.0keV_scaled.mat')
max(X)
min(X)
dlmwrite('ngc2300_MOS1_evt_0.5-8.0keV.txt', X, 'delimiter', '\t');






