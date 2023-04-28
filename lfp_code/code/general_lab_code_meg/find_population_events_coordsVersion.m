function popEvents = find_population_events_coordsVersion(allSpkTmsbyCell, coords, varargin)
% function popEvents = find_population_events(allSpkTmsbyCell, coords, speedOpt)
%
%   PURPOSE:
%        The purpose of this function is to identify population events that
%        would serve as candidate replay events. See Hwaun & Colgin 2019
%        and Pfeiffer & Foster 2013 for previous implementation of a
%        similar method.
%
%   INPUT:
%       allSpkTmsbyCell = cell(1,numUnits) that contains the spike times
%           for each unit in CA1 (generally want at least 20 for this type
%           of analysis)
%       coords = coords struct from read_in_coords
%       speedOpt = 1 to use times when speed < 5 cm/s (default), 2 for no 
%           speed check
%
%   OPTIONS:
%       Internally, there are several options that can be changed. These
%       include: how the binned firing rate is smoothed, the standard
%       deviations for detecting peaks and making the edge cuts, the number
%       of units that must be active in a bin for edge refining, the min
%       and max duration for population events, and the minimum number of
%       cells that must be active in each event.
%
%   OUTPUT:
%       popEvents = struct with times of the population events
%           (:,1) = start times
%           (:,2) = end times
%
%
% MM Donahue
% 04/2020
% Colgin Lab

%% OPTIONS

if nargin == 3
    speedOpt = varargin{1};
else
    speedOpt = 1;
end %check for input

% For getting/smoothing binned firing rate
binSz = 1/1000; %1 ms, as in Pfeiffer & Foster 2013

gWinStd = 10/1000; %10 ms, as in Pfeiffer & Foster 2013 & BN 2022
gWinDur = 50/1000; %in Berners-Lee et al. 2022, 100 ms
% gWinDur = gWinStd * 6;

% Detecting peaks and edges
pkCut = 3; %SD above the mean, as in Pfeiffer & Foster 2013 (and Hwaun & Colgin 2019)
edgeCut = 0; %as in Pfeiffer & Foster 2013 (and Hwaun & Colgin 2019)

% Min number of active units in first and last bin
minSpks = 1; %2 in Pfeiffer & Foster 2013

% Duration min and max
minDur = 50/1000; %50 ms
maxDur = 2000/1000; %2000 ms in Pfeiffer and Foster 2013

% Min number of cells that must be active in an event
% Suggested: 5 cells or 10% of all cells, whichever is higher
% minCell(1) = 5; %option 1
% totCell = length(allSpkTmsbyCell); %determine how many cells
% minCell(2) = floor(0.1 * totCell); %how much is 10% of them, option 2
% cellCrit = max(minCell); %see what the criteria should be
cellCrit = 5;


%% CHECK THAT THERE ARE AT LEAST 20 CELLS

if length(allSpkTmsbyCell) < 20
    warning('Less than 20 simultaneously recorded cells!')
end

%% GET SPIKE TIMES FROM ALL UNITS
fprintf('\t\t\t\tGetting spike times at run threshold\n')

if speedOpt == 1
    runThresh = 5; %cm/s -- threshold for movement
    instRs = get_runspeed(coords);
    smRs = smooth_runspeed(instRs);
    
    allSpkTms = []; %initalize
    usedSpkTmsbyCell = cell(1,length(allSpkTmsbyCell));
    
    for u = 1:length(allSpkTmsbyCell)
        tmpSpkTms = allSpkTmsbyCell{1,u};
        
        for st = 1:length(tmpSpkTms)
            rsInd = match(tmpSpkTms(st), smRs(:,1));
            if smRs(rsInd,2) < runThresh
                allSpkTms = [allSpkTms; tmpSpkTms(st)];
                usedSpkTmsbyCell{u} = [usedSpkTmsbyCell{u} tmpSpkTms(st)];
            end %run speed under thresh
        end %spike times
    end %units
    
elseif speedOpt == 2
    try
    allSpkTms = horzcat(allSpkTmsbyCell{:});
    catch
          allSpkTms = vertcat(allSpkTmsbyCell{:});
    end
    usedSpkTmsbyCell = allSpkTmsbyCell;
end %speed option

%% GET BINNED POPULATION FIRING RATE AND ZSCORE IT
fprintf('\t\t\t\tGetting binned population firing rate\n')

gWinStd = gWinStd / binSz; %convert based on bin size
gWinDur = gWinDur / binSz;
gKrnl = gausskernel(gWinDur, gWinStd);

startTm = coords(1,1);
endTm = coords(end,1);
edges = startTm:binSz:endTm;

spkCnts = zeros(1, length(edges));
if length(edges) < 2^16
    %can just use histcounts
    keyboard
else %cannot use histcounts
    for st = 1:length(allSpkTms)
        binInd = find(allSpkTms(st) >= edges, 1, 'Last'); 
        spkCnts(binInd) = spkCnts(binInd) + 1;
    end %all spike tms
end %whether or not we can use histcounts

binFr = spkCnts / binSz; %convert to Hz
smBinFr = conv(binFr, gKrnl, 'same'); %smooth it

zFr = zscore(smBinFr);

%% IDENTIFY POTENTIAL POPULATION EVENTS DURING IMMOBILE PERIODS
fprintf('\t\t\t\tIdentifying potential population events\n')

potEvInds = []; %initialize for storing potential events

findPeaks = find(zFr >= pkCut); %find where it peaks above 3
indDiffs = diff(findPeaks); %find the ind difference between the peaks...
findPeaks(indDiffs == 1) = []; %...and get rid of the ones that are right next to each other

for i = 1:length(findPeaks)
    pkInd = findPeaks(1,i);
    
    startEdge = find(zFr(1:pkInd) < edgeCut, 1, 'Last'); %find where it crossed edge thresh before...
    endEdge = find(zFr(pkInd:end) < edgeCut, 1, 'First');
    endEdge = endEdge + pkInd - 1; %...and after the peak
    
    if ~isempty(startEdge) && ~isempty(endEdge) %if the event was within this immobile period and not cut off
       potEvInds = [potEvInds; startEdge endEdge];
    end %not cut off
end %peaks

%% GET RID OF DETECTED EVENTS THAT SHARE EDGES
fprintf('\t\t\t\tRemoving redundant events\n')

i = 2;
while i <= size(potEvInds,1)
    if potEvInds(i,2) == potEvInds(i-1,2) %if these two are the same
        potEvInds(i,:) = []; %delete one
    else
        i = i + 1;
    end
end

%% REFINE EDGES SO THAT FIRST AND LAST ESTIMATION BINS CONTAIN A MIN # OF SPKS
fprintf('\t\t\t\tRefining event edges\n')

for i = 1:size(potEvInds,1)
    startInd = potEvInds(i,1);
    endInd = potEvInds(i,2);
    
    eventSpkCnts = spkCnts(startInd:endInd);
    
    rStartInd = find(eventSpkCnts, 1, 'First'); %refined start ind
    rStartInd = startInd + rStartInd - 1;
    rEndInd = find(eventSpkCnts, 1, 'Last'); %refined end ind
    rEndInd = startInd + rEndInd - 1;
    
    potEvInds(i,:) = [rStartInd rEndInd];
end %potential events

%% GET EVENT TIMES FROM INDS
fprintf('\t\t\t\tGetting event times\n')

potEvents = []; %master potential events array for storing times

for i = 1:size(potEvInds,1)
    startEdge = potEvInds(i,1);
    endEdge = potEvInds(i,2);
   
    startTm = edges(startEdge);
    endTm = edges(endEdge) + binSz;
    
    potEvents = [potEvents; startTm endTm];
end %i (potential events)

%% MAKE SURE EVENTS ARE 50-2000 ms
fprintf('\t\t\t\tChecking event durations\n')

badInds = []; %for catching inds for the events that aren't long enough

for i = 1:size(potEvents,1)
    dur = diff(potEvents(i,:));
    
    if dur < minDur || dur > maxDur
        badInds = [badInds; i]; %store ind for deletion
    end
end %inds

potEvents(badInds,:) = []; %delete bad ones

%% MAKE SURE AT LEAST 5 DIFFERENT CELLS OR 10% OF ALL CELLS ARE ACTIVE IN EVENT
fprintf('\t\t\t\tMaking sure enough cells are active in event\n')

badInds = []; %initialized

for i = 1:size(potEvents,1)
    actCell = 0; %initialize
    
    startTm = potEvents(i,1);
    endTm = potEvents(i,2);
    
    for u = 1:length(usedSpkTmsbyCell) %for each unit
        tmpSpkTms = usedSpkTmsbyCell{1,u}; %pull out the spike times
        spkInt = find(tmpSpkTms > startTm & tmpSpkTms < endTm); %find if this cell fired in the potential event
        
        if ~isempty(spkInt) %if it did
            actCell = actCell + 1; %count it
        end
    end %u
    
    if actCell < cellCrit %if after checking all the cells, not enough fired
        badInds = [badInds; i]; %store the bad inds
    end
end %i

potEvents(badInds,:) = []; %delete them

%% CREATE FINAL OUTPUT STRUCT WITH A BETTER SOUNDING NAME

popEvents = potEvents; %rename the struct for all the events that passed all criteria
% fprintf('DONE!\n')

end %function

