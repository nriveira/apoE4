function [coFirRm, coFirField] = get_2d_cofiring(u1SpkTms, u2SpkTms, coords)
% function [coFirRm, coFirField] = get_2d_cofiring(u1SpkTms, u2SpkTms, coords)
%
% PURPOSE:
%   The purpose of this function is to determine the cofiring ratemap and
%   field between two cells. The process is similar to O'Neill et al.,
%   2008. This function defines cofiring as firing within the same 100 ms
%   time bin. It determines the cofiring ratemap and cofiring field using
%   the same method as determining the individual cell's firing ratemaps
%   and place fields. Recommended to only use cell pairs that have a place
%   field similarity of at least 0.1.
%
% INPUT:
%   u1SpkTms = spike times from the first cell
%   u2SpkTms = spike times from the second cell
%   coords = coordinates of rat
%
% OUTPUT:
%   coFirRm = cofiring ratemap for the two cells
%   coFirField = cofiring field for the two cells, with subfields:
%       - coFirField.inds = x and y index values for field boundaries
%       - coFirField.cms = x and y cm values for field boundaries
%       - coFirField.pkFr = peak in-field cofiring rate (in Hz)
%       - coFirField.pkCoords = spatial index value for peak cofiring rate
%
% MMD
% 01/2022
% Colgin Lab

%% INITIALIZE

xBnds = [0 100];
yBnds = [0 100];
spatBinSz = 4; %cm
plotOrNot = 0;
velFilt = 1;
durCrit = 1;

minPkFr = 2; %Hz
minPfArea = 5; %cm^2

timeStep = 100/1000; %100 ms
timeBins = coords(1,1):timeStep:coords(end,1);

%% CALCULATE

if length(timeBins) > 2^16
    keyboard
end %more than max bins allowable by histcounts - shouldn't be but just in case!

u1BinSpks = histcounts(u1SpkTms, timeBins);
u2BinSpks = histcounts(u2SpkTms, timeBins);

spkProd = u1BinSpks .* u2BinSpks; %determine when they cofired - will be 0 if not, 1+ if yes
coFirEvs = timeBins(spkProd > 0); %time bins when there was co-firing

coFirRm = get_2d_ratemap(coFirEvs, coords, xBnds, yBnds, spatBinSz, plotOrNot, velFilt, durCrit); %ratemap method
coFirRm = smooth_2d_ratemap(coFirRm);
coFirField = get_2d_pfs(coFirRm, spatBinSz, minPkFr, minPfArea, plotOrNot);

% if ~isempty(coFirField)
%     fieldArea = polyarea(coFirField.inds(1,:), coFirField.inds(2,:)); %bin^2 area
%     if fieldArea/prod(size(coFirRm)) > 0.7
%         flag = 'Spatially unselective';
%     end %spatial selctive
% else
%     flag = 'Empty';
% end %whether there is a cofiring field to check

end %function