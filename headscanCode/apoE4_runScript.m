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

saveDir = '../figures/ADGrant/sample headscans_v2';
fileLoc = dir('../headscanData/iteration1');
for d = 3:length(fileLoc)
    plot_headscans([fileLoc(d).folder filesep fileLoc(d).name], param)
    figure(3); 
    saveas(gcf, [saveDir filesep extractBefore(fileLoc(d).name, 'DLC_resnet50')], 'png')
    saveas(gcf, [saveDir filesep extractBefore(fileLoc(d).name, 'DLC_resnet50')], 'eps')
    saveas(gcf, [saveDir filesep extractBefore(fileLoc(d).name, 'DLC_resnet50')], 'svg')
end

%% Plot all headscans
%       Calls plot_all_headscans function that plots all extracted data
tic
plot_all_headscans
toc