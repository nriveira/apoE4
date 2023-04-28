function popEvents = find_population_events(allSpkTmsbyCell, coords, speedOpt)
% function popEvents = find_population_events(allSpkTmsbyCell, coords)
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
%       speedOpt = 1 to use identified immobile periods
%           (find_immobile_periods), 2 to use times when speed < 5 cm/s, 3
%           for no speed check
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

% For getting/smoothing binned firing rate
binSz = 1/1000; %1 ms, as in Pfeiffer & Foster 2013

gWinStd = 10/1000; %10 ms, as in Pfeiffer & Foster 2013  & BN 2022
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
minCell(1) = 5; %option 1
totCell = length(allSpkTmsbyCell); %determine how many cells
minCell(2) = floor(0.1 * totCell); %how much is 10% of them, option 2
cellCrit = max(minCell); %see what the criteria should be

%% CHECK THAT THERE ARE AT LEAST 20 CELLS

if length(allSpkTmsbyCell) < 20
    warning('Less than 20 simultaneously recorded cells!')
end

%% IDENTIFY PERIODS OF IMMOBILITY/SPEED < 5 CM/S

% Pfeiffer and Foster 2013 used times when rat's velocity was <5cm/s
if speedOpt == 1
    %     fprintf('\tIdentifying periods of immobility\n')
    
    immTimes = find_immobile_periods(coords);
elseif speedOpt == 2
    %     fprintf('\tIdentifying speed < 5 cm/s\n')
    
    %same method as velocity filtering for ratemaps
    runThresh = 5; %cm/s -- threshold for movement
    
    instRs = get_runspeed(coords);
    [mwRs, mwInds] = mw_avg(instRs(:,2),15,7,1); % moving window average the runspeed to 0.5 s windows with 0.23 s steps
    mwRs = [instRs(1:mwInds(1)-1) mwRs instRs(mwInds(end)+1:length(instRs))]; %make the length of the moving window version match
    velFiltBnry = zeros(1,length(mwRs));
    velFiltBnry(mwRs<runThresh) = 1; %1 = less than 5 cm/s
    
    chunks = bwconncomp(velFiltBnry,4); %find where multiple consecutive frames of slow speed are
    
    immTimes = [];
    for ch = 1:chunks.NumObjects
        curChunk = chunks.PixelIdxList{ch};
        
        startInd = curChunk(1,1);
        endInd = curChunk(end,1);
        
        startTm = coords(startInd,1); %convert to seconds
        endTm = coords(endInd,1);
        
        immTimes = [immTimes; startTm endTm]; %save it for output
        
    end %chunks
elseif speedOpt == 3
    immTimes = [coords(1,1) coords(end,1)];
end %speed option

%% GET SPIKE TIMES FROM ALL UNITS

allSpkTms = []; %initalize

allSpkTms = horzcat(allSpkTmsbyCell{:});

% for u = 1:length(allSpkTmsbyCell)
%     tmpSpkTms = allSpkTmsbyCell{1,u};
% 
%     allSpkTms = [allSpkTms; tmpSpkTms];
% 
% end %units

%% GET BINNED POPULATION FIRING RATE
% fprintf('\tGetting binned population firing rate\n')

gWinStd = gWinStd / binSz; %convert based on bin size
gWinDur = gWinDur / binSz;
gKrnl = gausskernel(gWinDur, gWinStd);

binFrbyImm = cell(1,size(immTimes,1)); %by immobile period
binFrbyImmUS = cell(1,size(immTimes,1)); %unsmoothed

for it = 1:size(immTimes,1)
    startTm = immTimes(it,1);
    endTm = immTimes(it,2);
    
    totTime = endTm - startTm; %total time immobile in s
    if totTime == 0
        binFr = [];
    else
        edges = startTm:binSz:endTm;
        if length(edges) > 2^16
            tmpBinCnts = zeros(1, length(edges));
            
            for st = 1:length(allSpkTms)
                if allSpkTms(st) < startTm || allSpkTms(st) > endTm
                    continue
                end %if the spike time doesn't occur within this imm  period
                
                stInd = match(allSpkTms(st), edges);
                tmpBinCnts(stInd) =  tmpBinCnts(stInd) + 1; %add this to the count for this bin
            end %spike times
            
        else %there aren't too many bins
            tmpBinCnts = histcounts(allSpkTms, edges);
        end %too many edges
        
        binFr = tmpBinCnts * (1/binSz); %convert to Hz
    end %there is time
    
    binFrbyImmUS{1,it} = binFr;
    
    if totTime ~= 0
        smBinFr = conv(binFr, gKrnl, 'same'); %smooth it
        binFrbyImm{1,it} = smBinFr;
    else
        binFrbyImm{1,it} = [];
    end %is there any time
    
end %all immTimes

%% IDENTIFY MEAN AND STD FROM SMOOTHED FR BINS

allBinFr = []; %initialize
for it = 1:size(immTimes,1)
    allBinFr = [allBinFr binFrbyImm{1,it}]; %combine all of the bins
end %it

meanAll = mean(allBinFr); %get mean
SDAll = std(allBinFr); %get SD

%% IDENTIFY POTENTIAL POPULATION EVENTS DURING IMMOBILE PERIODS
% fprintf('\tIdentifying potential population events\n')

potEventsbyImm = cell(1,size(immTimes,1)); %1 = ind, 2 = times
for it = 1:size(immTimes,1)
    pullFr = binFrbyImm{1,it}; %pull out the bin firing rate for this event
    
    zFr = (pullFr - meanAll) ./ SDAll; %zscore using the overall mean and SD
    
    findPeaks = find(zFr > pkCut); %find where it peaks above 3
    indDiffs = diff(findPeaks); %find the ind difference between the peaks...
    findPeaks(indDiffs == 1) = []; %...and get rid of the ones that are right next to each other
    
    for i = 1:length(findPeaks)
        pkInd = findPeaks(1,i);
        
        startEdge = find(zFr(1:pkInd) < edgeCut, 1, 'Last'); %find where it crossed edge thresh before...
        endEdge = find(zFr(pkInd:end) < edgeCut, 1, 'First');
        endEdge = endEdge + pkInd - 1; %...and after the peak
        
        if ~isempty(startEdge) && ~isempty(endEdge) %if the event was within this immobile period and not cut off
            potEventsbyImm{1,it} = [potEventsbyImm{1,it}; startEdge endEdge]; %store inds
        end %not cut off
    end
end %immTimes

%% GET RID OF DETECTED EVENTS THAT SHARE EDGES
% fprintf('\tRemoving redundant events\n')

for it = 1:size(immTimes,1)
    i = 2;
    while i <= size(potEventsbyImm{1,it},1)
        if potEventsbyImm{1,it}(i,2) == potEventsbyImm{1,it}(i-1,2) %if these two are the same
            potEventsbyImm{1,it}(i,:) = []; %delete one
        else
            i = i+1;
        end
    end
end %it

%% REFINE EDGES SO THAT FIRST AND LAST ESTIMATION BINS CONTAIN A MIN # OF SPKS
% fprintf('\tRefining event edges\n')

for it = 1:size(immTimes,1)
    pullFr = binFrbyImmUS{1,it};
    for i = 1:size(potEventsbyImm{1,it},1)
        startInd = potEventsbyImm{1,it}(i,1);
        endInd = potEventsbyImm{1,it}(i,2);
        
        eventFr = pullFr(startInd:endInd);
        
        rStartInd = find(eventFr >= minSpks/binSz, 1, 'First'); %refined start ind
        rStartInd = startInd + rStartInd - 1;
        rEndInd = find(eventFr >= minSpks/binSz, 1, 'Last'); %refined end ind
        rEndInd = startInd + rEndInd - 1;
        
        potEventsbyImm{1,it}(i,:) = [rStartInd rEndInd];
        
    end %potential events
end %it

%% GET EVENT TIMES FROM INDS
% fprintf('\tGetting event times\n')

potEvents = []; %master potential events array
for it = 1:size(immTimes,1)
    immStartTm = immTimes(it,1);
    immEndTm = immTimes(it,2);
    
    edges = immStartTm:binSz:immEndTm;
    for i = 1:size(potEventsbyImm{1,it},1)
        startEdge = potEventsbyImm{1,it}(i,1);
        endEdge = potEventsbyImm{1,it}(i,2);
        %
        %         startTm = immTimes(it,1) + (startEdge * binSz) - binSz; %get it in seconds
        %         endTm = immTimes(it,1) + (endEdge * binSz) - binSz; %get it in seconds
        
        startTm = edges(startEdge) - binSz; %bin rounding error
        endTm = edges(endEdge) + binSz;
        
        potEvents = [potEvents; startTm endTm];
        
    end %i (potevential events)
end %immTimes


%% MAKE SURE EVENTS ARE 50-2000 ms
% fprintf('\tChecking event durations\n')

badInds = []; %for catching inds for the events that aren't long enough

for i = 1:size(potEvents,1)
    dur = potEvents(i,2) - potEvents(i,1);
    
    if dur < minDur || dur > maxDur
        badInds = [badInds; i]; %store ind for deletion
    end
    
end %inds

potEvents(badInds,:) = []; %delete bad ones

%% MAKE SURE AT LEAST 5 DIFFERENT CELLS OR 10% OF ALL CELLS ARE ACTIVE IN EVENT
% fprintf('\tMaking sure enough cells are active in event\n')

badInds = []; %initialized

for i = 1:size(potEvents,1)
    actCell = 0; %initialize
    
    startTm = potEvents(i,1);
    endTm = potEvents(i,2);
    
    for u = 1:length(allSpkTmsbyCell) %for each unit
        tmpSpkTms = allSpkTmsbyCell{1,u}; %pull out the spike times
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

