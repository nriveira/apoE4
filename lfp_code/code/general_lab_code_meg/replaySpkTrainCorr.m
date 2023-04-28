function corrSP = replaySpkTrainCorr(evSpkTms, uPkPos)
% function corrSP = replaySpkTrainCorr(evSpkTms, uPkPos)
% 
% PURPOSE:
%   Correlates the spike times within a replay event with the unit's peak
%   firing position on the track using circular-linear regression.
% 
% INPUT:
%    seqSpkTms = cell that is 1 x numCells in size, with each cell
%       containing all spike times from one unit during sequence
%           NOTE: Within each cell, spike times should be ordered as a
%           #spikes x 1 column. Otherwise, the correlation coefficient will
%           be incorrect.
%   uPkPos = 1 x numCells with the peaking firing position in degrees
% 
% OUTPUT:
%   corrSP = correlation between spike times and peak firing position.
% 
% MMD
% 11/2021
% Colgin Lab

if max(size(evSpkTms)) ~= length(uPkPos)
    error('Number of cells must be the same')
end

if max(uPkPos) > 2*pi
    uPkPos = deg2rad(uPkPos);
end %need to convert to rad

trackDiameter = 100; %cm
if nargin > 3
    trackDiameter = varagin{3};
end

spkTmsForCorr = []; %initialize
pkPosForCorr = [];

for u = 1:length(uPkPos)
        spkTmsForCorr = [spkTmsForCorr evSpkTms{u}'];
        pkPosForCorr = [pkPosForCorr repelem(uPkPos(u), length(evSpkTms{u}))];
end %units

if length(spkTmsForCorr) > 1
    corrSP = circ_corrcl(pkPosForCorr, spkTmsForCorr);
else
    corrSP = [];
end %there was anything to correlate

end %function