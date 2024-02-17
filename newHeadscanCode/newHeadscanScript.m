param.track_in_cm = 50;
param.fps = 30;
param.running_thresh = 6;
param.radial_thresh = 5;
param.pause_thresh = 6;
param.v_max = 50;      

saveDir = '../figures/new headscans/sample headscans';
fileLoc = dir('../headscanData/iteration1');
curDir = pwd;
lastInnerRadius = 300;
lastInnerCenter = [600, 350];
lastOuterRadius = 340;
lastOuterCenter = [600, 350];

%% Manually label track dimensions/coords without rat
for d = 3:length(fileLoc)
    % Look up the frame 
    lookupVal = split(extractBefore(fileLoc(d).name, 'DLC_resnet50'), '_');
    track_info(d-2).session = extractBefore(fileLoc(d).name, 'DLC_resnet50');
    videoDir = dir(['../headscanData/centerfit_videos' filesep lookupVal{1} filesep lookupVal{2}]);
    v = VideoReader([videoDir(3).folder filesep videoDir(3).name]);
    frame = read(v, 1);
    figure(1); clf; imshow(frame);
    inner = drawcircle('Color', 'b', 'FaceAlpha', 0.4, 'Center', lastInnerCenter, 'Radius', lastInnerRadius); pause;
    outer = drawcircle('Color', 'c', 'FaceAlpha', 0.4, 'Center', lastOuterCenter, 'Radius', lastOuterRadius); pause;

    % Assign all of values for next iteration
    lastInnerRadius = inner.Radius;
    lastInnerCenter = inner.Center;
    lastOuterRadius = outer.Radius;
    lastOuterCenter = outer.Center;
    track_info(d-2).radius = (outer.Radius + inner.Radius)./2;
    track_info(d-2).center = (inner.Center + outer.Center)./2;
    track_info(d-2).innerRadius = inner.Radius;
    track_info(d-2).innerCenter = inner.Center;
    track_info(d-2).outerRadius = outer.Radius;
    track_info(d-2).outerCenter = outer.Center;
end

%% Analyze points to find velocity, radial distance, and radial velocity
for d = 3:length(fileLoc)
    dlc_struct = get_dlc_struct_NHS([fileLoc(d).folder filesep fileLoc(d).name], track_info(d-2).center, track_info(d-2).radius);
    subplot(2,2,1); plot(dlc_struct(1).cartVelocity); title('Velocity (Cartesian)')
    subplot(2,2,2); plot(dlc_struct(1).predicted_r - 50); title('Radial Position')
    subplot(2,2,3); plot(dlc_struct(1).polar_coord); title('Angular Position')
    subplot(2,2,4); plot(dlc_struct(1).dadt); title('da/dt')
end