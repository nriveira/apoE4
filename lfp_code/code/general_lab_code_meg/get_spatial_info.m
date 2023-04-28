function spatInfo = get_spatial_info(rateMap, timePerBin)
% function spatInfo = get_spatial_info(rateMap, timePerBin)
%
% PURPOSE:
%   Computes spatial information. For formula reference, see Colgin et al.
%   2009, Mably et al. 2016, or the OG Skaggs et al. 1993/1996.
%
% INPUT:
%   rateMap = n x m matrix of firing rate by position bin
%   timePerBin = n x m matrix of time spent in each position bin
%
% OUTPUT:
%   spatInfo = computed spatial information for unit
%
% Code is adapted from code written by AM and CZ for the Colgin Lab.
%
% MM Donahue
% Colgin Lab
% 5/2021

if size(rateMap) ~= size(timePerBin)
   warning('Both inputs must be the same size!') 
   keyboard
end

rateMap = reshape(rateMap,1,[]); %reshape to linearize
timePerBin = reshape(timePerBin,1,[]); 

nanInds = find(isnan(timePerBin)); %find any bins that the rat did not visit
rateMap(nanInds) = [];
timePerBin(nanInds) = [];

nanInds = find(isnan(rateMap)); %find any bins that the rat did not visit
rateMap(nanInds) = [];
timePerBin(nanInds) = [];

numBins = length(rateMap);
meanFr = mean(rateMap); %same as lambda in formula
totTime = sum(timePerBin(:));

spatInfo = 0;
for i = 1:numBins
        
    pBin = timePerBin(i)/totTime; %prob of being in i bin
    lambI = rateMap(i); %mean firing rate of unit in i bin
   
    binSpatInfo = pBin * (lambI/meanFr) * log2(lambI/meanFr);
    
    if ~isnan(binSpatInfo)
         spatInfo = spatInfo + binSpatInfo; %in AM and CZ code, ignore where bin Fr = 0, time per bin = 0
    end %real number

end %bins


end %function