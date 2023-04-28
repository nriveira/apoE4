function [popEventInds, varargout] = find_population_events_spkRstr(spkRstr, Fs, varargin)
% function [popEventInds, popEventTms] = find_population_events_spkRstr(spkRstr, Fs, rstrTms)
%
%   PURPOSE:
%        The purpose of this function is to identify population events that
%        would serve as candidate replay events. See Hwaun & Colgin 2019
%        and Pfeiffer & Foster 2013 for previous implementation of a
%        similar method.
%
%   INPUT:
%       spkRstr = matrix that contains the time of each spike from each
%           cell in cell ID X sampling time point for each cell content
%               RECOMMENDED: only spikes from when the rat was traveling <5
%               cm/s should be included. This code does not velocity filter
%               internally.
%       Fs = sampling frequency of the spike raster
%               RECOMMENDED: use 1000 Hz sampling frequency to match the 
%               1 ms bin size used in most studies.
%       rstrTms = array of times that correspond to each raster time point
%               NOTE: only necessary if you are getting the popEventTms 
%               output  
%
%   OUTPUT:
%       popEventInds = inds (from spike raster) of the population events
%       popEventTms = tms of population events
%           (:,1) = start inds/times
%           (:,2) = end inds/times
%
% MM Donahue
% 04/2020
% Colgin Lab

%% CHECK

if nargout == 2 && nargin ~= 3
    error('Must include raster time mat to get population event times')
end

if nargin == 3
    rstrTms = varargin{1};
end

%% OPTIONS

binSz = 1/Fs;

gWinStd = 10/1000; %10 ms, as in Pfeiffer & Foster 2013 & BN 2022
gWinDur = 50/1000; %in Berners-Lee et al. 2022, 100 ms
% gWinDur = gWinStd * 6;

% Detecting peaks and edges
pkCut = 3; %SD above the mean, as in Pfeiffer & Foster 2013 (and Hwaun & Colgin 2019)
% edgeCut is 0 (see lines 86-101)

% Duration min and max
minDur = 50/1000; %50 ms
maxDur = 2000/1000; %2000 ms in Pfeiffer and Foster 2013

minBins = minDur * Fs;
maxBins = maxDur * Fs;

cellCrit = 5;

%% GET BINNED POPULATION FIRING RATE AND ZSCORE IT
% fprintf('\t\t\t\tGetting binned population firing rate\n')

gWinStd = gWinStd / binSz; %convert based on bin size
gWinDur = gWinDur / binSz;
gKrnl = gausskernel(gWinDur, gWinStd);

spkRstrFull = full(spkRstr);
popRstr = sum(spkRstrFull,1);
clear spkRstrFull %don't need it don't keep it

binFr = popRstr / binSz; %convert to Hz
smBinFr = conv(binFr, gKrnl, 'same'); %smooth it

zFr = zscore(smBinFr);

%% IDENTIFY POTENTIAL POPULATION EVENTS DURING IMMOBILE PERIODS
% fprintf('\t\t\t\tIdentifying potential population events\n')

potEvInds = []; %initialize for storing potential events

[~, pkLocs] = findpeaks(zFr, 'MinPeakHeight', pkCut); %find where it peaks above 3
zci = find(diff(sign(zFr))); %find where zscored firing rate crosses 0

for i = 1:length(pkLocs)
    pkInd = pkLocs(i);

    startEdgeInd = find(zci < pkInd, 1, 'Last');
    startEdge = zci(startEdgeInd);

    endEdgeInd = find(zci > pkInd, 1, 'First');
    endEdge = zci(endEdgeInd);
     
    if isempty(startEdge) || isempty(endEdge) 
      continue %to check next event 
    end %not cut off
     potEvInds = [potEvInds; startEdge endEdge];
end %peaks

%% GET RID OF DETECTED EVENTS THAT SHARE EDGES
% fprintf('\t\t\t\tRemoving redundant events\n')

i = 2;
while i <= size(potEvInds,1)
    if potEvInds(i,2) == potEvInds(i-1,2) %if these two are the same
        potEvInds(i,:) = []; %delete one
    else
        i = i + 1;
    end
end

%% REFINE EDGES SO THAT FIRST AND LAST ESTIMATION BINS CONTAIN A MIN # OF SPKS
% fprintf('\t\t\t\tRefining event edges\n')

for i = 1:size(potEvInds,1)
    startInd = potEvInds(i,1);
    endInd = potEvInds(i,2);
    
    eventRstr = popRstr(startInd:endInd); %population raster across this event
    
    rStartInd = find(eventRstr, 1, 'First'); %refined start ind
    rStartInd = startInd + rStartInd - 1;
    rEndInd = find(eventRstr, 1, 'Last'); %refined end ind
    rEndInd = startInd + rEndInd - 1;
    
    potEvInds(i,:) = [rStartInd rEndInd];
end %potential events

%% MAKE SURE EVENTS ARE 50-2000 ms
% fprintf('\t\t\t\tChecking event durations\n')

indDiffs = potEvInds(:,2) - potEvInds(:,1);
goodEvs = indDiffs>= minBins & indDiffs<=maxBins;
potEvInds = potEvInds(goodEvs,:);

%% MAKE SURE AT LEAST 5 DIFFERENT CELLS ARE ACTIVE IN EVENT
% fprintf('\t\t\t\tMaking sure enough cells are active in event\n')

goodInds = []; %initialized

for i = 1:size(potEvInds,1)
    startInd = potEvInds(i,1);
    endInd = potEvInds(i,2);
    
    evSpkRstr = full(spkRstr(:,startInd:endInd));
    actCell = nnz(sum(evSpkRstr,2));

    if actCell >= cellCrit
         goodInds = [goodInds; i]; %store the good inds
    end
end %i

popEventInds = potEvInds(goodInds,:); 

%% GET TIMES IF NEEDED

if nargout == 2
    popEventTms = zeros(size(popEventInds));
    for i = 1:size(popEventInds,1)
        popEventTms(i,:) = [rstrTms(popEventInds(i,1)) rstrTms(popEventInds(i,2))];
    end %i - events
end %get times bc they have been requested

varargout{1} = popEventTms;

end %function

