%%  Plot headscans
%       Calls analyze_headscans and plots figures for individual DLC
%       results
param.track_in_cm = 50;
param.bodypart = 1;
param.median_window = 15;
param.fps = 30;
param.running_thresh = 6;
param.radial_thresh = 3;
param.pause_thresh = 6;
param.v_max = 25;                                                                                                                                                                                                                                                          

plot_headscans("../headscanData/iteration1/rat387_d1s2DLC_resnet50_headscansFeb24shuffle1_100000.csv", param)
title('Sample Headscan detection')

%% Plot all headscans
%       Calls plot_all_headscans function that plots all extracted data
tic
plot_all_headscans
toc