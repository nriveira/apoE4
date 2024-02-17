%%  Plot headscans
%       Calls analyze_headscans and plots figures for individual DLC
%       results
param.track_in_cm = 50;
param.median_window = 15;
param.fps = 30;
param.running_thresh = 6;
param.radial_thresh = 5;
param.pause_thresh = 6;
param.v_max = 50;      

saveDir = '../figures/poster';
fileLoc = dir('../headscanData/iteration1');
for d = 3:length(fileLoc)
    extractBefore(fileLoc(d).name, 'DLC_resnet50');
    plot_headscans([fileLoc(d).folder filesep fileLoc(d).name], param)
    figure(3);
    saveas(gcf, [saveDir filesep extractBefore(fileLoc(d).name, 'DLC_resnet50')], 'png')
    saveas(gcf, [saveDir filesep extractBefore(fileLoc(d).name, 'DLC_resnet50')], 'epsc')
end

%% Plot all headscans
%       Calls plot_all_headscans function that plots all extracted data
tic
plot_all_headscans
toc