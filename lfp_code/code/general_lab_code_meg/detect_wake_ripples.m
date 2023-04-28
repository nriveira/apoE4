function ripInfo = detect_wake_ripples(allLfpStruct, coords, speedCheck)
% function ripInfo = detect_wake_ripples(allLfpStruct, coords)
%
%   PURPOSE:
%        The purpose of this function is to detect sharp wave ripples while
%        the animal is awake. This function combines ripples detected
%        across all tetrodes in CA1. For reference, see Karlsson & Frank
%        2009.
% 
%   INPUT:
%       allLfpStruct = cell(1,numCA1tets), where each cell contains the
%           lfpStruct from all tetrodes that were in CA1. Each lfpStruct
%           should contain the fields:
%               - allLfpStruct{1,x}.data*
%               - allLfpStruct{1,x}.ts*
%               - allLfpStruct{1,x}.Fs*
%               - allLfpStruct{1,x}.fn*
%               - allLfpStruct{1,x}.deltaFiltLfp**
%               - allLfpStruct{1,x}.thetaFiltLfp**
%               - allLfpStruct{1,x}.ripFiltLfp**
%           *read_in_lfp
%           **filter_lfp_for_rip_detection
%       coords = coords struct from read_in_coords
%       speedCheck = option for how to check to animal's speed during
%           detected events
%               1 = identify immobile periods (see find_immobile_periods),
%               recommended for open field/big square/free exploration data
%               2 = check run speed is under maximum, recommended for
%               circle track SWR identification while eating rewards etc.
%               0 = don't check either one (not recommended - I don't know
%               why you would do this)
%
%   OPTIONS:
%       Internally, there are some options that can be changed, such as the
%       min and max duration for events to be included and the time window
%       around each event to evaluate the theta/delta ratio.
%
%   OUTPUT:
%       ripInfo.inds = inds that define each ripple
%       ripInfo.tms = times that define each ripple
%           *for both, (:,1) = onset, (:,2) = offset
%
%
% MM Donahue
% 04/2020
% Colgin Lab

%% OPTIONS
durThresh = [0.05 0.25]; 

winLength = 5; %s, for evaluating theta/delta ratio (and speed if using speed check 2)

pkCut = 3; %this is flexible. Karlsson, Jadhav use 3, Cheng used 3 OR 6 depending on what they called the event
pkTime = 0.015; 
edgeCut = 0; 

sampFreq = 2000; %lfp sampling frequency

%% CHECK LFP STRUCT

for c = 1:length(allLfpStruct)
    if length(allLfpStruct{1,c}.data) ~= length(allLfpStruct{1,1}.data)
       fprintf('\tCheck CSC %d\n', c)
       keyboard
    end
end %tet for check


%% INITIALIZE THINGS

c = 1; %for setting up ripBinary
ripBinary = zeros(1,length(allLfpStruct{1,c}.data));

pkTimeInd = pkTime * sampFreq; %convert peak time based on sampling frequecy

gaussWin = 0.032; %in second; 32 ms for smoothing, as in Jadhav et al. 
gaussSd = 0.004; %in second; 4 ms for smoothing, as in Karlsson & Frank 2009
gaussSd = gaussSd * sampFreq; %convert based on sampling freqeuncy
kernel = gausskernel(gaussWin,gaussSd);

%% IDENTIFY PEAKS

for c = 1:length(allLfpStruct)
    lfpStruct = allLfpStruct{1,c};
    hilbLfp = abs(hilbert(lfpStruct.ripFiltLfp)); %Hilbert transform
    smLfp = conv(hilbLfp, kernel, 'same'); %
%     ripPow = zscore(abs(hilbert(lfpStruct.ripFiltLfp)));

    ripPow = zscore(smLfp);

    findPeaks = find(ripPow >= pkCut);
    
    indDiffs = [0; diff(findPeaks)];
    allPeaks = findPeaks; %keep for next step
    findPeaks(indDiffs == 1) = []; %delete the ones that are right next to each other
    
    %% CHECK IF THEY PEAK FOR MIN TIME
    
    badInds = []; %initialize
    
    for i = 1:length(findPeaks)
        pkInd = findPeaks(i);
        endPeakAboveCut = find(ripPow(pkInd:end) < pkCut, 1, 'First'); %how many inds away does it dip below our peak cut-off
        
        if endPeakAboveCut < pkTimeInd
           badInds = [badInds; i]; %save for deletion
        end
    end %peaks to check time
    
    findPeaks(badInds) = [];
    
    %% GET EDGES AROUND EACH PEAK
    
    ripEdges = []; %initialize
    
    for i = 1:length(findPeaks)
        pkInd = findPeaks(i,1);
        
        startEdge = find(ripPow(1:pkInd) < edgeCut, 1, 'Last');
        endEdge = find(ripPow(pkInd:end) < edgeCut, 1, 'First');
        endEdge = endEdge + pkInd - 1;
        
        if ~isempty(startEdge) & ~isempty(endEdge) %SWRS are not cut off
            ripEdges(i,:) = [startEdge endEdge]; 
        end
        
    end %all found peaks
    
    %% GET RID OF MULTIPLES
    
    i = 2;
    while i <= size(ripEdges,1)
        if ripEdges(i,2) == ripEdges(i-1,2)
            ripEdges(i,:) = [];
        else
            i = i+1;
        end
    end
    
    %% GET RIP TIMES
    try
     ripEdgeTimes = lfpStruct.ts(ripEdges);
    catch
        if ripEdges(1,1) == 0
            ripEdges(1,:) = [];
            ripEdgeTimes = lfpStruct.ts(ripEdges);
        end
    end
    
    if size(ripEdgeTimes,2) ~= 2
        ripEdgeTimes = ripEdgeTimes';
    end

    %% CHECK THETA/DELTA RATIO DURING SWR
    
    deltaAmpTs = abs(hilbert(lfpStruct.deltaFiltLfp));
    thetaAmpTs = abs(hilbert(lfpStruct.thetaFiltLfp));
    
    thetaDeltaRatio = thetaAmpTs ./ deltaAmpTs;
    
    thetaDeltaRatio(thetaDeltaRatio == Inf) =  max(thetaDeltaRatio(thetaDeltaRatio~=Inf)); %get rid of infinites
    thetaDeltaRatio(thetaDeltaRatio == 0) = 0.0001; %get rid of zeros
    
    avg = nanmean(thetaDeltaRatio);
    stdev = nanstd(thetaDeltaRatio);
    
    thetaDeltaRatio = (thetaDeltaRatio - avg) ./ stdev;
    
    ripCtrInds = round(mean(ripEdges,2));
    
    i = 1;
    while i <= size(ripCtrInds)
        
        try
            windowTDR = mean(thetaDeltaRatio(ripCtrInds(i) - round(lfpStruct.Fs * (winLength/2)):ripCtrInds(i) + round(lfpStruct.Fs * (winLength/2))));
        catch %if the window start is before the recording started
            if ripCtrInds(i) - round(lfpStruct.Fs * (winLength/2)) < 1
                windowTDR = mean(thetaDeltaRatio(1:ripCtrInds(i) + round(lfpStruct.Fs * (winLength/2)))); %just go as early as you can
            elseif ripCtrInds(i) + round(lfpStruct.Fs * (winLength/2)) > length(thetaDeltaRatio) %if the window end is later than the recording
                windowTDR = mean(thetaDeltaRatio(ripCtrInds(i) - round(lfpStruct.Fs * (winLength/2)):end)); %just go as late as you can
            end
        end
        if windowTDR > 0
            ripCtrInds(i) = [];
            ripEdgeTimes(i,:) = [];
            ripEdges(i,:) = [];
        else
            i = i+1;
        end
    end
    
    %% ADD IT TO RIPPLE BINARY
    
    for i = 1:size(ripEdges,1)
        ripBinary(ripEdges(i,1):ripEdges(i,2)) = 1;
    end %i

end %for all CSCs

%% PULL RIPPLES FROM THE RIP BINARY

chunks = bwconncomp(ripBinary,4); %find connect components

ripInds = [];
ripTms = [];

for SWR = 1:chunks.NumObjects
    startInd = chunks.PixelIdxList{1,SWR}(1);
    endInd = chunks.PixelIdxList{1,SWR}(end);
    try
    startTm = lfpStruct.ts(startInd);
    catch
        keyboard
    end
    try
    endTm = lfpStruct.ts(endInd);
    catch
        keyboard
    end
    
    ripInds = [ripInds; startInd endInd];
    ripTms = [ripTms; startTm endTm];
end %all detected SWRs

%% COMBINE EVENTS WITHIN 25 MS OF EACH OTHER

i = 2;
while i < size(ripInds,1)
    if (ripTms(i,1) - ripTms(i-1,2)) <= .025 %from Ernie
        ripInds(i-1,2) = ripInds(i,2);
        ripTms(i-1,2) = ripTms(i,2);
        ripInds(i,:) = [];
        ripTms(i,:) = [];
    else
        i = i+1;
    end
end

%% CHECK DURATION
    
badInds = []; %initialize
for i = 1:size(ripInds,1)
    dur = ripTms(i,2) - ripTms(i,1);
    
    if dur < durThresh(1) | dur > durThresh(2)
        badInds = [badInds; i];
    end %if
end %i

ripInds(badInds,:) = [];
ripTms(badInds,:) = [];

%% SPEED CHECK - VERIFY THAT THE ANIMAL WAS IMMOBILE DURING SWR
    
if speedCheck == 1
    immTimes = find_immobile_periods(coords);
    
    badInds = []; %initialize
    
    for i = 1:size(ripInds,1)
        startTm = ripTms(i,1);
        endTm = ripTms(i,2);
        immCheck = 0; %initialize it as bad
        for it = 1:size(immTimes,1)
            if startTm > immTimes(it,1) & startTm < immTimes(it,2) & endTm > immTimes(it,1) & endTm < immTimes(it,2)
                immCheck = 1; %it occured during an immobile period so it's good
                break
            end %if
            
        end %immobile times
        
        if immCheck == 0
            badInds = [badInds; i]; %save for deletion
        end
        
    end %inds
    
    ripInds(badInds,:) = [];
    ripTms(badInds,:) = [];
    
%% SPEED CHECK - VERIFY ANIMAL WAS MOVING <5 CM/S
    
elseif speedCheck == 2
    
    maxSpd = 5; %cm/s
    instRunSpd = get_runspeed(coords);
    
    badInds = []; %initialize
    
    for i = 1:size(ripInds,1)
        startTm = ripTms(i,1);
        
        evWin = [startTm - (winLength / 2) startTm + (winLength /2)]; %time window to evaluate
%         spdWin = [find(instRunSpd(:,1) >= evWin(1), 1, 'Last') find(instRunSpd(:,1) <= evWin(2), 1, 'First')
        spdWin = [find(instRunSpd(:,1) >= evWin(1)& instRunSpd(:,1) <= evWin(2))];
        evSpd = mean(instRunSpd(spdWin,2));
        
        if evSpd > maxSpd
            badInds = [badInds; i];
        end %too fast
        
    end %events (i)
    
    ripInds(badInds,:) = [];
    ripTms(badInds,:) = [];
    
end %which speedcheck

%% SET UP OUTPUT

ripInfo.tms = ripTms;
ripInfo.inds = ripInds;

end %fxn