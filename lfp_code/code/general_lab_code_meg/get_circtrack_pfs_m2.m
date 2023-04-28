function pf = get_circtrack_pfs_backwardsExpansion(rateMap, spatBinSz, varargin)
% function pf = get_circtrack_pfs_backwardsExpansion(rateMap, spatBinSz, trackDiameter)
% 
% PURPOSE:
%   This function uses the Bieri et al. 2014 method to detect place fields
%   (modified for the circle track). This det

pf = []; %initialize

rmBinSz = 360/length(rateMap); %ratemap bin size in deg
degBinCtrs = rmBinSz/2:rmBinSz:360;

minPfLen = 18; %degrees, apprx. 15 cm in Bieri et al. 2014
contigBins = round(minPfLen/rmBinSz);

minFr = 2.5; %Hz,average firing rate across ratemap had to exceed

if max(rateMap) < minFr
    return
end

stdMinFr = 1; %Hz, standard minimum firing rate
pkFrCell = max(rateMap);

minFrUse = min([stdMinFr 0.1*pkFrCell]); %firing rate above the lower of either 1 Hz or 10% of the peak firing rate of the cell

abvPeak = rateMap >= minFrUse;
CC = bwconncomp(abvPeak, 4); %get connected components

inds = [];
for p = 1:CC.NumObjects
    if length(CC.PixelIdxList{1,p}) >= contigBins
        tmpInds = [CC.PixelIdxList{1,p}(1) CC.PixelIdxList{1,p}(end)];
        inds = [inds; tmpInds];
    end %if long enough
end %potential place fields

%check if there's anything that wraps around 0 deg
if CC.PixelIdxList{1,1}(1) == 1 && CC.PixelIdxList{1,end}(end) == length(rateMap)
    tmpAllInds = [CC.PixelIdxList{1,end}' CC.PixelIdxList{1,1}'];
    if length(tmpAllInds) >= contigBins
        %delete old potential pfs to replace with mroe inclusive one
        oldPfs = [];
        tmpOld = [oldPfs find(inds(:,1) == tmpAllInds(1))];
        if ~isempty(tmpOld)
            oldPfs = [oldPfs tmpOld];
        end
        tmpOld =  find(inds(:,2) == tmpAllInds(end));
        if ~isempty(tmpOld)
            oldPfs = [oldPfs tmpOld];
        end
        inds(oldPfs,:) = [];
        
        inds = [inds; tmpAllInds(1) tmpAllInds(end)];
    end %long enough
end %potential pf around 0

newInds = []; %initialize
delInds = [];
for p = 1:size(inds,1)
    if inds(p,2) > inds(p,1)
        tmpInds = inds(p,1):inds(p,2);
    else
        tmpInds = [inds(p,1):length(rateMap) 1:inds(p,2)];
    end
    pullRawRm = rateMap(tmpInds);
    badInds = find(pullRawRm==0);
    if isempty(badInds)
        continue
    end %whether there is anything to process
    
    delInds = [delInds; p];
    tmpPotInds = tmpInds(1:badInds(1)-1);
    if length(tmpPotInds) > contigBins
        newInds = [newInds;  tmpPotInds(1) tmpPotInds(end)];
    end
    for i = 1:length(badInds)-1
        tmpPotInds = tmpInds(badInds(i)+1:badInds(i+1)-1);
        if length(tmpPotInds) > contigBins
            newInds = [newInds;  tmpPotInds(1) tmpPotInds(end)];
        end
    end %inds
    tmpPotInds = tmpInds(badInds(end)+1:end);
    if length(tmpPotInds) >= contigBins
        newInds = [newInds;  tmpPotInds(1) tmpPotInds(end)];
    end
end %potenital pfs

inds(delInds,:) = [];
inds = [inds; newInds];

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