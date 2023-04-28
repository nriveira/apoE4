function dir = detect_sequence_events(pxn, spkRstr, Fs, radPos, radBinCtrs, bothDir)
% function [seqInds, seqTms, seqSlopes] = detect_sequence_events(pxn, spkRstr, Fs, radPos, radBinCtrs, bothDir)
%
% PURPOSE:
%   To detect sequence events using the method stated in Zheng et al. 2021.
%       Works for circle track data, including free runs where the rat
%       crosses 0 degrees (unlike in DMS task). For code reference, see
%       DetectSequenceEvents_cz.
%
% INPUTS:
%   pxn = output from BayesianDecoder function
%   spkRstr = SAME spike raster used as input for Bayesian Decoder function
%       to get pxn
%   Fs = sampling frequency of the spike raster
%   radPos = radial position of the rat over the same time period as pxn
%   radBinCtrs = radial position of each bin for pxn
%   bothDir = 1 for both forward and reverse sequences, 0 for forward only
%
% OUTPUTS:
%	dir.dirInd = forward (1) or reverse (2) event
%      .seq:
%           inds = inds from the Bayesian decoder for each seq
%           tms = times for each seq
%           slope = slope, in rad/s
%           pxn = pxn across seq
%           com = center of mass of distribution
%           calphase = calculated phase of the regression line
%           actPos = actual position of the rat
%
% OPTIONS:
%   CHECK FUNCTION FOR THRESHOLDS! Set to match those set in Zheng et al.
%       2021.
%
% MMD
% 7/2021
% Colgin Lab
keyboard
%% OPTIONS/INITIALIZE

maxJumpThr = 1.4; %rad - "estimated positions between adjacent time bins did not exceed 1.4 rad"
timeWin = 6; %6 subsequent time bins with spike
timeStep = 1; %bin step
distanceThr = 0.07; %"distance between first and last estimated position within seqeunce was more or equal to 0.07 rad"

bayesWin = 40/1000; %40 ms
bayesStep = 10/1000; %10 ms
sampRate = 20000; %Hz - spike sampling rate for spike raster input to Bayesian decoder

minSpk = 5; %5 spikes in an event
minCell = 3; %3 cells must spike in event

ppToFitThr = 0.35; %rad, "at least 60% of  the total post prob needed to be no mroe than 0.35 rad away from the fitted traj line"
propFitThr = 0.6;

fitToActThr = 0.35; %rad, minimal distance between the fitted trajectory and rats' true postion had to be less than or equal to 0.35 rad

if max(radPos(:,2)) > 2*pi
    radPos(:,2) = deg2rad(radPos(:,2));
end

if bothDir == 1
    seqTimes = cell(2,1);
    seqSlopes = cell(2,1);
    dirs = [1 2];
else
    seqTimes = [];
    seqSlopes = [];
    dirs = 1;
end %both dir or not

dirNames = {'Forward', 'Reverse'};

onset = cell(2,1);
offset = cell(2,1);

%% IDENTIFY POTENTIAL EVENTS FROM PXN

for ii = 1:timeStep:size(pxn,2)
    range = ii:ii+timeWin-1;
    if max(range)>size(pxn,2)
        break
    end %min range fits in pxn
    
    emptyBin = isnan(pxn(1,range));
    if any(emptyBin)
        continue %to next ind
    end %any nan
    
    %     [~, decodedPos] = max(pxn(:,range)); %decoded position based on max likely
    
    com = NaN(1,length(range)); %initialize
    for t = 1:length(range)
        com(t) = circ_mean(radBinCtrs', pxn(:,range(t))); %weighted circular mean based on prob distribution
    end %time bins
    
    jump = nan(1,length(com)-1);
    for ij = 1:length(jump)
        jump(ij) = circ_dist(com(ij+1), com(ij));
    end %get jump
    
    maxJump = nanmax(abs(jump));
    
    seqDist = circ_dist(com(end), com(1));
    distFromStart = circ_dist(com, com(1));
    
    if bothDir == 1
        if maxJump <= maxJumpThr && abs(seqDist) >= distanceThr
            if seqDist > 0
                dirInd = 1;
            else
                dirInd = 2;
            end %direction of seq - potential forward or rev?
            
            if dirInd == 1 && min(distFromStart) >= 0 || dirInd == 2 && max(distFromStart) <=0
                onset{dirInd} = cat(1,onset{dirInd},ii);
                offset{dirInd} = cat(1,offset{dirInd},range(end));
            end
        end %qualifies
    else %forward only
        dirInd = 1;
        if maxJump <= maxJumpThr && seqDist >= distanceThr && min(distFromStart) >= 0
            %                 seq = [seq; onOff];
            onset{dirInd} = cat(1,onset{dirInd},ii);
            offset{dirInd} = cat(1,offset{dirInd},range(end));
        end %meets crit
    end %both directions
end %pxn bins

%% COMBINE OVERLAPPING EVENTS

for dirInd = dirs %not combining forward events with reverse events
    if length(onset{dirInd}) > 1
        isi = onset{dirInd}(2:end)-offset{dirInd}(1:end-1);
        merge = find(isi <= 0);
        onset{dirInd}(merge+1) = [];
        offset{dirInd}(merge) = [];
    end
end %dir ind

%% CHECK CELLS/SPIKES

%make nsameple a odd number - done in BayesianDecoder function
if mod(Fs,2) == 0
    nsample = bayesWin*Fs+1;
else
    nsample = bayesWin*Fs;
end

tRange = floor(nsample/2)+1:bayesStep*Fs:size(spkRstr,2)-floor(nsample/2); %time range
for dirInd = dirs
    badInds = []; %initialize
    
    for s = 1:length(onset{dirInd})
        rangeOn = tRange(onset{dirInd}(s))-floor(nsample/2); %range from pxn
        try
            rangeOff = tRange(offset{dirInd}(s))+floor(nsample/2)-1;
        catch
            if offset{dirInd}(s) > length(tRange)
                rangeOff = size(spkRstr,2);
            end
        end %time range doesn't fit
        
        pullRstr = full(spkRstr(:,rangeOn:rangeOff));
        
        numSpksByCell = sum(pullRstr,2);
        
        if sum(pullRstr(:)) < minSpk || length(numSpksByCell(numSpksByCell~=0)) < minCell %if not good
            badInds = [badInds; s];
        end %meet spike and cell crit
    end %all seq
    
    if ~isempty(onset{dirInd})
        onset{dirInd}(badInds) = [];
        offset{dirInd}(badInds) = [];
    end %anything to delete
end %dirInd

%% GET TIMES

onsetTm = cell(2,1);
offsetTm = cell(2,1);

[nWin, winStartInds] = find_num_windows(size(spkRstr,2), bayesWin*sampRate, bayesStep*sampRate);
if winStartInds(end)+bayesWin*sampRate < size(spkRstr,2)

    nWin = nWin + 1;
    winStartInds(end+1) = winStartInds(end)+bayesStep*sampRate; %okay
end
winStartTms = radPos(1,1) + winStartInds/sampRate - 1/sampRate;
winEndTms = winStartTms + bayesWin;

for dirInd = dirs
    badInds = []; %initialize
    
    for s = 1:length(onset{dirInd})
        if offset{dirInd}(s) > nWin
            badInds = [badInds s];
        else
            tmpOn = onset{dirInd}(s);
            tmpOff = offset{dirInd}(s);
            
            onTm = winStartTms(tmpOn);
            offTm = winEndTms(tmpOff);
            onsetTm{dirInd} = [onsetTm{dirInd} onTm];
            offsetTm{dirInd} = [offsetTm{dirInd} offTm];
        end %within actual nWin
    end %seq
    
    onset{dirInd}(badInds) = [];
    offset{dirInd}(badInds) = [];
end %dir ind

% %% CHEK RUN SPEED
% 
% runSpd = get_runspeed(coords);
% smRunSpd = smooth_runspeed(runSpd);
% 
% for dirInd = dirs
%     badInds = []; %initialize
%     for s = 1:length(onset{dirInd})
%         onTm = onsetTm{dirInd}(s);
%         offTm = offsetTm{dirInd}(s);
%         
%         startInd = find(smRunSpd(:,1) <= onTm, 1, 'Last');
%         endInd = find(smRunSpd(:,1) <= offTm, 1, 'Last');
%         
%         seqSpd = smRunSpd(startInd:endInd, 2);
%         if min(seqSpd) < runThresh %if speed drops below threshold during the sequence
%             badInds = [badInds s];
%         end %speed check
%     end %sequences
%     
%     onset{dirInd}(badInds) = [];
%     onsetTm{dirInd}(badInds) = [];
%     offset{dirInd}(badInds) = [];
%     offsetTm{dirInd}(badInds) = [];
%     
% end %dirInd

%% CHECK SUFFICIENT % PXN IS CLOSE TO FITTED LINE / CHECK FITTED CLOSE ENOUGH TO TRUE POS
keyboard
seqSlopes = cell(2,1); %initialize
for dirInd = dirs
    badInds = []; %re-initialize
    
    for s = 1:length(onset{dirInd})
        
        tmpOn = onset{dirInd}(s);
        tmpOff = offset{dirInd}(s);
        
        pullPxn = pxn(:,tmpOn:tmpOff);
        
        com = NaN(1,size(pullPxn,2)); %initialize
        for t = 1:length(com)
            com(t) = circ_mean(radBinCtrs', pullPxn(:,t)); %weighted circular mean based on prob distribution
        end %time bins
        com = wrapTo2Pi(com);
        
        maxVals = max(pullPxn);
        if sum(maxVals==1/size(spkRstr,1)) > 0
            keyboard
        end
        timeVals = 0:bayesStep:(size(pullPxn,2)*bayesStep)-bayesStep;
        
        %         com(maxVals==1/size(spkRstr,1)) = []; %no better than chance
        %         timeVals(maxVals==1/size(spkRstr,1)) = [];
        
        beta = CircularRegression(timeVals, com);
        slope = beta(1);
        seqSlopes{dirInd} = [seqSlopes{dirInd}; slope];
        
        calphase = beta(1)*timeVals + beta(2);
        
        if dirInd == 1 && slope < 0
            badInds = [badInds s];
            continue
        elseif dirInd == 2 && slope > 0
            badInds = [badInds s];
            continue
        end %slope facing right way
        
        sumPxn = 0;
        for t = 1:length(timeVals)
            fitPos = calphase(t); %fit position
            
            fitPosMin = fitPos - ppToFitThr;
            fitPosMax = fitPos + ppToFitThr;
            
            [~,minInd] = min(abs(circ_dist(radBinCtrs, fitPosMin)-0));
            [~,maxInd] = min(abs(circ_dist(radBinCtrs, fitPosMax)-0));
            
            if maxInd > minInd
                sumPxn = sumPxn + sum(pullPxn(minInd:maxInd,t));
            else
                sumPxn = sumPxn + sum(pullPxn(minInd:end,t)) + sum(pullPxn(1:maxInd,t)); %fit crosses 0
            end  
        end %time vals
        
        sumAll = sum(pullPxn(:));
        propWithDist = sumPxn/sumAll; %proportion of pxn within distance
        
        if propWithDist < propFitThr
            badInds = [badInds s];
            continue
        end %not high enough %
        
        %     keyboard
        %     figure;
        %     plot(timeVals, radDecPos)
        % %     y_reg = 2*pi*para(1)*timeVals + para(2);
        %     hold on;
        %     plot(timeVals, calphase)
        %     ylim([0 2*pi])
        
        onTm = onsetTm{dirInd}(s); %time seq starts in s
        offTm = offsetTm{dirInd}(s);
        
        startInd = match(onTm, radPos(:,1));
        endInd = match(offTm, radPos(:,1));
        keyboard 
        %circ_mean
        actPos = mean(radPos(startInd:endInd,2));
        
        minDist = min(abs(circ_dist(calphase, actPos)));
        
        if minDist > fitToActThr
            badInds = [badInds s];
        end %doesn't meet min dist crit
    end %seq
    
    onset{dirInd}(badInds) = [];
    onsetTm{dirInd}(badInds) = [];
    offset{dirInd}(badInds) = [];
    offsetTm{dirInd}(badInds) = [];
    seqSlopes{dirInd}(badInds) = [];
end %dirInd

%% PREPARE OUTPUT

for dirInd = dirs
    seqCntr = 0;
    dir(dirInd).name = dirNames{dirInd};
    for i = 1:length(onset{dirInd})
        seqCntr = seqCntr + 1;
        
        dir(dirInd).seq(seqCntr).dirInd = dirInd;
        dir(dirInd).seq(seqCntr).inds = [onset{dirInd}(i) offset{dirInd}(i)];
        dir(dirInd).seq(seqCntr).tms = [onsetTm{dirInd}(i) offsetTm{dirInd}(i)];
        dir(dirInd).seq(seqCntr).pxn = pxn(:,onset{dirInd}(i):offset{dirInd}(i));
        dir(dirInd).seq(seqCntr).slope = seqSlopes{dirInd}(i);
        
        %just going to re-do all these here, it's just easier on the brain
        
        startInd = match(dir(dirInd).seq(seqCntr).tms(1), radPos(:,1));
        endInd = match(dir(dirInd).seq(seqCntr).tms(2), radPos(:,1));
        
        actPos = mean(radPos(startInd:endInd,2));
        
        dir(dirInd).seq(seqCntr).actPos = actPos;
        
        com = NaN(1,size( dir(dirInd).seq(seqCntr).pxn,2)); %initialize
        for t = 1:length(com)
            com(t) = circ_mean(radBinCtrs',  dir(dirInd).seq(seqCntr).pxn(:,t)); %weighted circular mean based on prob distribution
        end %time bins
        com = wrapTo2Pi(com);
        
        dir(dirInd).seq(seqCntr).com = com;
        
        timeVals = 0:bayesStep:(length(com)*bayesStep)-bayesStep;
        beta = CircularRegression(timeVals, com);
        calphase = beta(1)*timeVals + beta(2);
        if beta(1) ~=  dir(dirInd).seq(seqCntr).slope
            keyboard
        end %check
        
        dir(dirInd).seq(seqCntr).calphase = calphase;
        
    end %seq
end %direction


end %function