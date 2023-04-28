function pf = get_circtrack_pfs_backwardExpansion(rateMap, varargin)
% function pf = get_circtrack_pfs_backwardExpansion(rateMap, minMeanFr, minPkFr, trackDiameter)
% 
% PURPOSE:
%   This function uses the Bieri et al. 2014 method to detect place fields
%   (modified for the circle track). 
% 
% INPUTS:
%   rateMap = rateMap across the circle track for this place field,
%       smoothed (in Bieri et al., smoothing w/ Gaussian kernel std = 25
%       cm)
%   minMeanFr = minimum mean firing rate across the ratemap to be
%       considered a potential place cell (in Bieri et al., 2.5 Hz).
%       Default = 0 Hz (not used).
%   minPkFr = minimum peak firing rate across the ratemap to be considered
%       a potential place cell (in Zheng et al. 2021, 1 Hz). Default = 1
%       Hz.
%   trackDiameter (optional) = diameter of circle track. Default = 100 cm
% 
% OUTPUT:
%   pf = struct for each place field with fields:
%       inds = ratemap indicies for field
%       radPos = position in degress of place field
%       pkFr = peak firing rate in field
%       pkPos = position in degress of peak firing in field
% 
% MMD
% Colgin Lab
% 06/2022

pf = []; %initialize

rmBinSz = 360/length(rateMap); %ratemap bin size in deg
degBinCtrs = rmBinSz/2:rmBinSz:360; %centers of spatial bins, in deg

minMeanFr = 0;
if nargin > 1
    minMeanFr = varargin{1};
end %check for track diameter

if mean(rateMap) < minMeanFr
    return
end %doesn't meet minimum

minPkFr = 1; %Hz
if nargin > 2
    minPkFr = varargin{2};
end %check for track diameter

if max(rateMap) < minPkFr
    return
end %doesn't meet minimum

trackDiameter = 100; %cm
if nargin > 3
    trackdiameter = varargin{3};
end %check for track diameter

minPfLen = 15 * (360/(trackDiameter*pi)); %degrees, apprx. 15 cm in Bieri et al. 2014
contigBins = round(minPfLen/rmBinSz); %needed to reach ~15 cm

stdMinFr = 1; %Hz, standard minimum firing rate
pkFrCell = max(rateMap);

minFrUse = min([stdMinFr 0.1*pkFrCell]); %firing rate above the lower of either 1 Hz or 10% of the peak firing rate of the cell

abvPeak = rateMap >= minFrUse;
CC = bwconncomp(abvPeak, 4); %get connected components

inds = [];
for p = 1:CC.NumObjects
    if length(CC.PixelIdxList{1,p}) >= contigBins
        tmpInds = [CC.PixelIdxList{1,p}(1)-1 CC.PixelIdxList{1,p}(end)+1];
        inds = [inds; tmpInds];
    end %if long enough
end %potential place fields

%check if there's anything that wraps around 0 deg
if CC.PixelIdxList{1,1}(1) == 1 && CC.PixelIdxList{1,end}(end) == length(rateMap)
    tmpAllInds = [CC.PixelIdxList{1,end}' CC.PixelIdxList{1,1}'];
    if length(tmpAllInds) >= contigBins
        %delete old potential pfs to replace with more inclusive one
        oldPfs = [];
        tmpOld = [oldPfs find(inds(:,1) == tmpAllInds(1)-1)];
        if ~isempty(tmpOld)
            oldPfs = [oldPfs tmpOld];
        end
        tmpOld =  find(inds(:,2) == tmpAllInds(end)+1);
        if ~isempty(tmpOld)
            oldPfs = [oldPfs tmpOld];
        end
        inds(oldPfs,:) = [];
        
        inds = [inds; tmpAllInds(1)-1 tmpAllInds(end)+1];
    end %long enough
end %potential pf around 0

inds(inds==0) = length(rateMap);
inds(inds==length(rateMap)+1) = 1;

for p = 1:size(inds,1)
    if inds(p,2) > inds(p,1)
        pf(p).inds = inds(p,1):inds(p,2);
        pf(p).radPos = degBinCtrs(pf(p).inds);
        
        [pf(p).pkFr, maxInd] = max(rateMap(pf(p).inds));
        pf(p).pkPos = degBinCtrs(pf(p).inds(1) + maxInd - 1);
    else
        pf(p).inds = [inds(p,1):length(rateMap) 1:inds(p,2)];
        pf(p).radPos = degBinCtrs(pf(p).inds);
        [pf(p).pkFr, maxInd] = max(rateMap(pf(p).inds));
        pf(p).pkPos = pf(p).radPos(maxInd);
    end %whether cross 0
end %place fields

end %function