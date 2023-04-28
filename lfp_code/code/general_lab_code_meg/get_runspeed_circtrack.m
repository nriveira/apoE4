function instRunSpd = get_runspeed_circtrack(radPos, varagin)
% function instRunSpd = get_runspeed_circtrack(radPos, trackDiameter)
%
% PURPOSE:
%   To calculate the run speed in cm/s from the radial position of the rat
%   on the circle track.
%
% INPUT:
%     radPos = n x 2 matrix where column 1 equals time stamps for each frame,
%              column 2 radial position of the rat (deg or rad acceptable)
%   trackDiameter (optional) = default is 100 cm (standard circular track
%       for Colgin Lab. Include if different track size is used.
%
% OUTPUT:
%   instRunSpd = n x 2 matrix where column 1 is frametimes and column 2 is instantaneous run speed for each frame
%
%
% MMD - edited from original get_runspeed code by JBT
% 09/2021
% Colgin Lab

trackDiameter = 100;
if nargin > 1
    trackDiameter = varagin{2};
end
degCmConv = trackDiameter*pi / 360;

if max(radPos(:,2)) < 2*pi
    radPos(:,2) = rad2deg(radPos(:,2)); %convert to deg
end

tBetFrames = mean(diff(radPos(:,1))); %find the average time between frames
frameRate = 1/tBetFrames;

% find degrees rat travelled between each frame
degDiff = abs(diff(radPos(:,2))); %difference in x coords between each frame and the frame 'frameOffset' beyond it
crossZeroDiff = find(degDiff > 300); %find candidate times when difference likely indicates rat cross 0-360 deg boundary
for i = 1:length(crossZeroDiff)
%     degDiff(crossZeroDiff(i)) = abs(rad2deg(circ_dist(deg2rad(radPos(crossZeroDiff(i),2)),deg2rad(radPos(crossZeroDiff(i)+1,2)))));
    degDiff(crossZeroDiff(i)) = 360 - degDiff(crossZeroDiff(i)); %actually the same thing as above
end %i

tmpRunSpds = degDiff * degCmConv; %convert to cm/frame
tmpRunSpds = tmpRunSpds .* frameRate; %convert to cm/s

instRunSpd(:,2) = [tmpRunSpds(1); tmpRunSpds]; %since run speed is the difference between frames, there is one less value here than there are for coords.
%                               To make things easy on ourselves, we'll assume the rat is moving the same speed from frame 0-1 as he is from 1-2.

instRunSpd(:,1) = radPos(:,1);


end %function