function corrSP = seqSpikeTrainCorr(seqSpkTms, uPkPos, actPos, vargin)
% function corrSP = seqSpikeTrainCorr(allSpkTms, uPkPos, actPos, trackDiamter)
% 
% PURPOSE:
%   Calculate the correlation coefficient between spike time and place
%       cell peak firing positions within 50 cm of the rat on the circle
%       track. See Feng, Silva, & Foster 2015 as reference and for more
%       detailed information.
% 
% INPUTS:
%   seqSpkTms = cell that is 1 x numCells in size, with each cell
%       containing all spike times from one unit during sequence
%           NOTE: Within each cell, spike times should be ordered as a
%           #spikes x 1 column. Otherwise, the correlation coefficient will
%           be incorrect.
%   uPkPos = 1 x numCells with the peaking firing position in degrees
%   trackDiameter (optional) = default is 100 cm (standard circular track
%       for Colgin Lab). Include if different track size is used in cm.
% 
% OUTPUT:
%   corrSP = correlation between spike times and peak firing position.
% 
% MMD
% 9/2021
% Colgin Lab

if max(size(seqSpkTms)) ~= length(uPkPos)
    error('Number of cells must be the same')
end

trackDiameter = 100; %cm
if nargin > 3
    trackDiameter = varagin{3};
end

degCmConv = trackDiameter*pi / 360;

normPos = 180; %normalized as if rat is at 180
shiftDeg = normPos - actPos;
shiftUPkPos = wrapTo360(uPkPos + shiftDeg);

spkTmsForCorr = []; %initialize
pkPosForCorr = [];

for u = 1:length(uPkPos)
    if shiftUPkPos(u) > (180 - 50/degCmConv) && shiftUPkPos(u) < (180 + 50/degCmConv) %within 50 cm of actual (now 180)
        spkTmsForCorr = [spkTmsForCorr seqSpkTms{u}'];
        pkPosForCorr = [pkPosForCorr repelem(shiftUPkPos(u), length(seqSpkTms{u}))];
    end %occurs within 50 cm
end %units

if length(spkTmsForCorr) > 1
    cMat = corrcoef(spkTmsForCorr, pkPosForCorr); %get correlation coefficient
    corrSP = cMat(2);
else
    corrSP = [];
end %there was anything to correlate


end %function