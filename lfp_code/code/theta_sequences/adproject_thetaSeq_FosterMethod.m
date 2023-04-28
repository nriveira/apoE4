%% EXAMPLE PIPELINE CODE FOR JAYANTH

% NOTE:
%   This is how I usually take my code and is how my FXS code is
%   formulated. There are some differences with how CZ and EH did it,
%   although the actual outcome is the same. I will try to make a note
%   about those differences.
%addpath 'E:\circular_track_data' 

cd('E:\Rat381\2022-10-24') %go into day folder

TTList = Readtextfile('TTList.txt'); %read in list of cells

close all;

clear coords coords_diff

%[radPos, coords] = read_in_circtrack_coords('VT1.nvt');

cd('E:\Rat381\2022-10-24\begin3')   
radPos = rat.day(1).begin(3).radPos;
coords = rat.day(1).begin(3).coords;
Ts = rat.day(1).begin(3).radPos(:,1); %convert and re-name


% this is John's code, adapted from Ernie's code, that reads in the circle
% track coordinates in both radial position (0-360deg) and x/y coordinates.
% For loading in W-track coords, I would recommend using the read_in_coords
% function. This will give you the x and y coordinates of the rat for each
% time point. You will then have to linearize it when you create the
% ratemap for each cell.

% CZ at some point normalized all of her data so that the start box was
% always at 0 deg/rad. You can check the function doall_data_pos2ang.m for
% how she did that.

% For the actual ratemap, I'm taking a lot of the code/verbiage from John's
% code, but repeating it here for clarity. You can look up his
% get_ratemap_circtrack if you'd like as well.

spatBinSz = 4; %spatial bin size of the ratemap, in degrees
spatBinEdges = 0:spatBinSz:360; %edges of our spatial bins
degBinCtrs = spatBinSz/2:spatBinSz:360; %centers of all spatial bins of the ratemap, in deg
radBinCtrs = deg2rad(degBinCtrs); %same, in radians

runThresh = 5; %cm/s, min speed rat must be running for a spike to be included in the ratemap
% durCrit = 150/1000; %150 ms; rat must spend a min amount of time in a bin for it to be included.
% This will not matter too much for linearized environments, as the rat
% will run through all bins multiple times, but if you are using a 2D
% ratemap, it will be important so you don't get erroneously high firing
% rates in some bins, so I just wanted to mention it.

gKrnl = gausskernel(5,2); %this was the kernel Ernie used for smoothing the ratemaps.
% This kernel uses a 10 bin window ([-5:5]) = 10*spatBinSz = 40 deg
% This kernel uses a 2 bin standard deviation = 2*spatBinSz = 8 deg
rmOffset = 0.0001; %Bayesian decoder gets messed up when there are ratemap bins = 0

timelimit = 120;   % from CZ: set the highest speed, unit is cm/s
runSpeed = speed2D(coords(:,2), coords(:,3), coords(:,1)); %velocity in cm/s
runSpeed(runSpeed>=timelimit) = 0.5*(runSpeed(circshift((runSpeed>=timelimit),-3)) + runSpeed(circshift((runSpeed>=timelimit),3)));
% you can also linearize the run speed if you'd like, CZ did this by
% calculating the angular velocity. See ratemap_AllLaps_v2.m (among others)
% also see get_runspeed function from John

slowInds = runSpeed<runThresh;
for i = 1:length(radPos)
    if radPos(i,2) < 0
        radPos(i,2) = radPos(i,2) + 360;
    end
end
filtRadPos = radPos(~slowInds,:);

%%
tpb = zeros(1,length(degBinCtrs)); %time per bin, denominator of ratemap
for sb = 1:length(spatBinEdges)-1
    binStart = spatBinEdges(sb);
    % this code uses bin edges, but you could do distance to bin centers as
    % well by doing something like [~, minInd] = min(circ_dist(radBinCtrs', deg2rad(filtRadPos(timeBin,2))));
    % more difficult for circular data since you need the circular distance
    binEnd = spatBinEdges(sb+1);
    frameInds = find(filtRadPos(:,2)>=binStart & filtRadPos(:,2)<binEnd);
    % +1 extra frame to go from 1/2 a frame before time-start to 1/2 a frame after time-end
    tpb(sb) = (length(frameInds)+1) * mean(diff(coords(:,1))); %frames x smapling rate
end %bins

%%
rateMaps = zeros(length(TTList), length(radBinCtrs)); %initialize for storing all ratemaps for all cells
uIDs = zeros(length(TTList),2); %useful for keep track of cells active in each event, not necessary
spksPerCellFilt = cell(length(TTList),1); %for decoding later - velocity filtered for replay
spksPerCell = cell(length(TTList),1); %not velocity filtered

for u = 1:length(TTList)
    [tetNum, clustNum] = get_unit_ID(TTList{u}); %again, optional but helpful
    uIDs(u,:) = [tetNum, clustNum];

    spkTms = Readtfile(TTList{u}); %get spike times from the .t files
    spkTms = spkTms ./ 10^4; %convert it from samples to seconds

    spkCntsByFrame = zeros(1,size(radPos,1));
    spdBySpk = zeros(length(spkTms),1); %initialize
    for st = 1:length(spkTms)
        tmpTm = spkTms(st);

        if tmpTm < radPos(1,1)
            continue;
        end

        frameInd = find(radPos(:,1)<=tmpTm, 1, 'Last');

        spkCntsByFrame(frameInd) = spkCntsByFrame(frameInd) + 1; %can be >1 because of sampling rate difs
        spdBySpk(st) = runSpeed(frameInd);
    end %spike times

    spkCntsByFrame = spkCntsByFrame(~slowInds); %velocity filter

    for sb = 1:length(spatBinEdges)-1
        binStart = spatBinEdges(sb); %this code uses bin edges, but you could do distance to bin centers as well
        binEnd = spatBinEdges(sb+1);
        frameInds = find(filtRadPos(:,2)>=binStart & filtRadPos(:,2)<binEnd);
        % +1 extra frame to go from 1/2 a frame before time-start to 1/2 a frame after time-end
        spkCnts(sb) = sum(spkCntsByFrame(frameInds));
    end %bins

    tmpMap = spkCnts ./ tpb; %calculate raw ratemap for this cell
    rateMaps(u,:) = conv_cir(tmpMap, gKrnl); %store the smoothed ratemap

    spksPerCell{u} = spkTms; %all spike times, not velocity filtered
    velFiltSpkTms = spkTms(spdBySpk<runThresh); %spikes under 5 cm/s for decoding
    spksPerCellFilt{u} = velFiltSpkTms;

end %unit


rateMaps = rateMaps + rmOffset; %can't use 0s in Bayesian Decoder
% NOW WE HAVE THE RATEMAPS FOR ALL THE PLACE CELLS COLLECTED FROM MULTIPLE
% TETRODES. WE ESSENTIALLY HAVE P(N/X) FOR THE BAYESIAN DECODER.


%% THETA SEQUENCES EXTRACTION
% WE WILL JUST INITIALIZE ALL THE VARIABLES REQUIRED HERE FOR GETTING THE THETA SEQUENCES
% Now we will detect theta sequences events from the decoded probability
% distribution. For this portion, I'm going to limit it to individual laps
% to make it easier.

% trackFileName = dir('*tracking*'); %where CZ stored info about lap times
% load(trackFileName.name, 'Ts_sample') %load in the times for the sample laps

Fs = 2000; %whatever seems appropriate, this is what CZ and EH used

bayesWin = 40/1000; %40 ms for decoding
bayesStep = 10/1000; %10 ms

% I'm going to use the wording from some code I wrote to detect events, but
% reproduce it here for clarity and keep it all together. It's called
% "detect_sequence_events" on the server. CZ's code does not allow for
% events that cross the 0-360 deg mark, since the rest box was there in the
% adjusted position for all of her experiments. In my FXS data, I needed
% there to be a possibility for the rat to cross that mark, since there was
% no rest box.

maxJumpThr = 1.4; %rad - "estimated positions between adjacent time bins did not exceed 1.4 rad"
timeWin = 6; %6 subsequent time bins with spike (based on bayes win and step) = MEANING 240 MS?
timeStep = 1; %bin step
distanceThr = 0.07; %"distance between first and last estimated position within seqeunce was more or equal to 0.07 rad"

runThresh = 5; %cm/s
% CZ did not institute a run speed when detecting events, so you could set
% this to 0 to eliminate this check

minSpk = 5; %5 spikes in an event
minCell = 3; %3 cells must spike in event

ppToFitThr = 0.35; %rad, "at least 60% of  the total post prob needed to be no mroe than 0.35 rad away from the fitted traj line"
propFitThr = 0.6;
fitToActThr = 0.35; %rad, minimal distance between the fitted trajectory and rats' true postion had to be less than or equal to 0.35 rad

forSeq = []; %initialize for storing all sequences from all sample laps - forward events
revSeq = []; %reverse events

spkRstr = make_spike_raster(spksPerCell, [Ts(1,1) Ts(length(Ts),1)], Fs); %time from leaving rest box to reaching reward
pxn = BayesianDecoder(spkRstr, rateMaps, bayesWin, bayesStep, Fs); %get pxn

    % ~~~~~~ First, do our initial detection of events with sequential
    % structure
theta_onset = cell(2,1); %two cells for forward and reverse
theta_offset = cell(2,1);
    for ii = 1:timeStep:size(pxn,2)
        range = ii:ii+timeWin-1;
        if max(range)>size(pxn,2)
            break
        end %fits within the pxn

        emptyBin = isnan(pxn(1,range));
        if any(emptyBin)
            continue %to next ind
        end %any nan

        %could use max prob as decoded bin, or com
        com = NaN(1,length(range)); %initialize
        for t = 1:length(range)
            com(t) = circ_mean(radBinCtrs', pxn(:,range(t))); %weighted circular mean based on prob distribution
        end %time bins
        
        %GET THE JUMP DISTANCE
        jump = nan(1,length(com)-1);
        for ij = 1:length(jump)
            jump(ij) = circ_dist(com(ij+1), com(ij));
        end %get jump

        maxJump = nanmax(abs(jump));

        seqDist = circ_dist(com(end), com(1)); %path length of sequence
        distFromStart = circ_dist(com, com(1)); %distance each bin from start of sequence

        %we can consider both forward and reverse theta sequences
        if maxJump <= maxJumpThr && abs(seqDist) >= distanceThr
            if seqDist > 0
                dirInd = 1; %forward events
            else
                dirInd = 2; %reverse event
            end %direction of seq - potential forward or rev?

            if dirInd == 1 && min(distFromStart) >= 0 || dirInd == 2 && max(distFromStart) <=0
                theta_onset{dirInd} = cat(1,theta_onset{dirInd},ii);
                theta_offset{dirInd} = cat(1,theta_offset{dirInd},range(end));
            end % doesn't swipe back for forward or forward for reverse
        end %qualifies
    end %time bin step

    %     ~~~~ Combine events that overlap
    for dirInd = 1:2 %not combining forward events with reverse events
        if length(theta_onset{dirInd}) > 1
            isi = theta_onset{dirInd}(2:end)-theta_offset{dirInd}(1:end-1);
            merge = find(isi <= 0);
            theta_onset{dirInd}(merge+1) = [];
            theta_offset{dirInd}(merge) = [];
        end %more than one event - check
    end %dir ind

    %      ~~~~~~ Check active cells and spikes - I'm going to delete the ones
    %      that don't meet the minimums, but you could save and store this info
    %      and only include events that meet mins later (will be easier to see
    %      how changing these parameters will change results)

    if mod(Fs,2) == 0
        nsample = bayesWin*Fs+1;
    else
        nsample = bayesWin*Fs;
    end %make nsample a odd number - done in BayesianDecoder function

    tRange = floor(nsample/2)+1:bayesStep*Fs:size(spkRstr,2)-floor(nsample/2); %time range for decoding

    for dirInd = 1:2
        badInds = []; %initialize for deletion
        for s = 1:length(theta_onset{dirInd})
            rangeOn = tRange(theta_onset{dirInd}(s))-floor(nsample/2); %range from pxn -> spk raster
            try
                rangeOff = tRange(theta_offset{dirInd}(s))+floor(nsample/2)-1;
            catch
                if theta_offset{dirInd}(s) > length(tRange)
                    rangeOff = size(spkRstr,2);
                end
            end %time range doesn't fit

            pullRstr = full(spkRstr(:,rangeOn:rangeOff)); %pull from master spike raster

            numSpksByCell = sum(pullRstr,2);

            if sum(pullRstr(:)) < minSpk || length(numSpksByCell(numSpksByCell~=0)) < minCell %if not good
                badInds = [badInds; s];
            end %meet spike and cell crit
        end %all seq

        if ~isempty(theta_onset{dirInd})
            theta_onset{dirInd}(badInds) = [];
            theta_offset{dirInd}(badInds) = [];
        end %anything to delete
    end %for/rev

    %     ~~~~~ Get the times of the events
    onsetTm = cell(2,1);
    offsetTm = cell(2,1);

    [nWin, winStartInds] = find_num_windows(size(spkRstr,2), bayesWin*Fs, bayesStep*Fs);
    %based on Fs and spike raster, how long would pxn be for only FULL 40ms windows
    if winStartInds(end)+bayesWin*Fs < size(spkRstr,2)
        nWin = nWin + 1;
        winStartInds(end+1) = winStartInds(end)+bayesStep*Fs; %okay
    end
    winStartTms = radPos(1,1) + winStartInds/Fs - 1/Fs; %equivalent times for pxn windows
    winEndTms = winStartTms + bayesWin;

    for dirInd = 1:2
        badInds = []; %initialize

        for s = 1:length(theta_onset{dirInd})
            if theta_offset{dirInd}(s) > nWin
                badInds = [badInds s];
            else
                tmpOn = theta_onset{dirInd}(s);
                tmpOff = theta_offset{dirInd}(s);

                onTm = winStartTms(tmpOn);
                offTm = winEndTms(tmpOff);
                onsetTm{dirInd} = [onsetTm{dirInd} onTm];
                offsetTm{dirInd} = [offsetTm{dirInd} offTm];
            end %within actual nWin
        end %seq

        theta_onset{dirInd}(badInds) = [];
        theta_offset{dirInd}(badInds) = [];
    end %dir ind

    % ~~~~~~ Check sufficient % of pxn is close to fitted line/check fitted
    % line is close enough to true pos/rat is running at least 5 cm/s

    seqSlopes = cell(2,1); %initialize
    seqCal = cell(2,1);
    for dirInd = 1:2
        badInds = []; %re-initialize

        for s = 1:length(theta_onset{dirInd})
            tmpOn = theta_onset{dirInd}(s);
            tmpOff = theta_offset{dirInd}(s);

            pullPxn = pxn(:,tmpOn:tmpOff);

            com = NaN(1,size(pullPxn,2)); %initialize
            for t = 1:length(com)
                com(t) = circ_mean(radBinCtrs', pullPxn(:,t)); %weighted circular mean based on prob distribution
            end %time bins
            com = wrapTo2Pi(com);

            timeVals = 0:bayesStep:(size(pullPxn,2)*bayesStep)-bayesStep; %for calculating regression
            bins2use = find(~isnan(sum(pullPxn)));

            [~, calphase, ~, ~, slope] = Cir_reg(pullPxn, radBinCtrs', timeVals, bins2use);
            seqSlopes{dirInd} = [seqSlopes{dirInd}; slope];
            seqCal{dirInd} = [seqCal{dirInd}; {calphase}];

            if dirInd == 1 && slope < 0 || dirInd == 2 && slope > 0
                badInds = [badInds s];
                continue
            end %slope facing right way

            sumPxn = 0; %to sum tp pxn close to fit line
            for t = 1:length(timeVals)
                fitPos = calphase(t); %fit position

                fitPosMin = fitPos - ppToFitThr; %loc in rad
                fitPosMax = fitPos + ppToFitThr;

                [~,minInd] = min(abs(circ_dist(radBinCtrs, fitPosMin)-0));
                [~,maxInd] = min(abs(circ_dist(radBinCtrs, fitPosMax)-0));

                if maxInd > minInd
                    sumPxn = sumPxn + sum(pullPxn(minInd:maxInd,t));
                else %crosses 0
                    sumPxn = sumPxn + sum(pullPxn(minInd:end,t)) + sum(pullPxn(1:maxInd,t)); %fit crosses 0
                end %whether or not crosses 0
            end %time vals

            sumAll = sum(pullPxn(:));
            propWithDist = sumPxn/sumAll; %proportion of pxn within distance

            if propWithDist < propFitThr
                badInds = [badInds s];
                continue
            end %not high enough %

            onTm = onsetTm{dirInd}(s); %time seq starts in s
            offTm = offsetTm{dirInd}(s);

            startInd = match(onTm, radPos(:,1)); %equivalent bin for this time for radPos
            endInd = match(offTm, radPos(:,1));
            actPos = circ_mean(deg2rad(radPos(startInd:endInd,2))); %mean of position

            minDist = min(abs(circ_dist(calphase,  actPos)));

            if minDist > fitToActThr
                badInds = [badInds s];
                continue
            end %doesn't meet min dist crit

            seqSpeed = runSpeed(startInd:endInd);
            if min(seqSpeed) < runThresh
                badInds = [badInds s];
            end %doesn't meet speed criteria

        end %s - seq

%         onset{dirInd}(badInds) = [];
%         onsetTm{dirInd}(badInds) = [];
%         offset{dirInd}(badInds) = [];
%         offsetTm{dirInd}(badInds) = [];
%         seqSlopes{dirInd}(badInds) = [];
%         seqCal{dirInd}(badInds) = [];
    end %dirInd

    %     ~~~~~~Prepare output
    lpForSeq = []; %initailize
    dirInd = 1;
    for s = 1:length(theta_onset{dirInd})
        lpForSeq(s).inds = [theta_onset{dirInd}(s) theta_offset{dirInd}(s)];
        %this is the inds of the pxn
        lpForSeq(s).tms = [onsetTm{dirInd}(s) offsetTm{dirInd}(s)];
        %time is absolute for whole day - doesn't matter if detected from
        %lap or day pxn
        lpForSeq(s).pxn = pxn(:,theta_onset{dirInd}(s):theta_offset{dirInd}(s));
        lpForSeq(s).slope = seqSlopes{dirInd}(s);
        lpForSeq(s).calphase = seqCal{dirInd}{s}';
    end %seq
    forSeq = cat(2, forSeq, lpForSeq); 
%     Depending on what you want to do or what makes more sense to you, you
%     may format the output differently. I usually store my data in a
%     struct organized rat - day - task (pre-run/trials/post-test) - lap
%     number - lap phase (sample/test), but EH and CZ did one struct with 
%     all events and included a field for rat number, day, lap, etc. You 
%     may also want to include information like whether it was correct or
%     incorrect, how far away the rat was if it was incorrect, etc.

    lpRevSeq = []; %initialize
    dirInd = 2;
    for s = 1:length(theta_onset{dirInd})
        lpRevSeq(s).inds = [theta_onset{dirInd}(s) theta_offset{dirInd}(s)];
        lpRevSeq(s).tms = [onsetTm{dirInd}(s) offsetTm{dirInd}(s)];
        lpRevSeq(s).pxn = pxn(:,theta_onset{dirInd}(s):theta_offset{dirInd}(s));
        lpRevSeq(s).slope = seqSlopes{dirInd}(s);
        lpRevSeq(s).calphase = seqCal{dirInd}{s}';
    end %seq
    revSeq = cat(2, revSeq, lpRevSeq);


%%
% Plot the theta sequence heatmap with x-axis = normalized time bins,
% y-axis = before and after the rat's current position

for i = 1:length(lpForSeq)
%     xtmp(i) = length(lpForSeq(i).pxn(1,:));
%     pxn_norm(i).pxn = normalize(lpForSeq(i).pxn, 'range');
    %ymid(i) = round((find(radPos(:,1)<=lpForSeq(i).tms(2), 1, 'Last') - find(radPos(:,1)<=lpForSeq(i).tms(1), 1, 'Last'))/2 + find(radPos(:,1)<=lpForSeq(i).tms(1), 1, 'Last'));
    ymid(i) = round((find(radPos(:,1)<=lpForSeq(i).tms(2), 1, 'Last') - find(radPos(:,1)>=lpForSeq(i).tms(1), 1, 'First'))/2 + find(radPos(:,1)>=lpForSeq(i).tms(1), 1, 'First'));
end
% xtmp = max(xtmp);

normTaxis = 0:0.05:1;

for i = 1:length(lpForSeq)
    %actPos = radBinCtrs(ymid(i));
    actPos_new = round(radPos(ymid(i),2));
    
    shiftVal = 45 - round(actPos_new/4);
    
%     if actPos_new > 0
%         shiftVal = 45 - round(actPos_new/4);
%     else
%         actPos_new = actPos_new + 360;
%         shiftVal = 45 - round(actPos_new/4);
%     end
    
    actPos_all(i).pos = actPos_new;
    pxn_ind = lpForSeq(i).pxn;
    in = ~isnan(sum(pxn_ind,1));
    pxn_ind(:,~in) = 1/length(radBinCtrs);
    txAx = 1:size(pxn_ind,2);
    txAx = txAx/max(txAx);
    normPxn = VFR(pxn_ind, txAx', normTaxis, size(pxn_ind,1), mean(diff(normTaxis)));
    %normPxn_all(i).pxn = flip(normPxn,1);
    normPxn_all(i).pxn = normPxn;
    shiftPxn(i).pxn = circshift(normPxn, shiftVal);%shift so rat's current position is at center
    shiftPxn(i).pxn = flip(shiftPxn(i).pxn,1);
end
% 
mean_pxn = mean(cat(4,shiftPxn(:).pxn),4);
imagesc(mean_pxn(36:54,:));
colorbar;
xticklabels(1/(length(xticks)):1/(length(xticks)+1):1);
yticklabels(10:-2.5:-10)
ylabel('10 position bins (40 degrees) before and after the rat''s current position');
xlabel('Normalized time between 0 and 1')
yline(10, 'k', 'Linewidth',2);


%%

posInds = 1:1:length(radPos);
tmInds = 1:size(pxn,2);
winMids = Ts(1,1)+((tmInds-1)*bayesStep)+(.5*bayesWin);

% Get rat's position for every bayes window
posMat = radPos(posInds,:);
allTms = winMids;
knownPts = posMat(:,1);
knownPos = posMat(:,2);
allPos = interp1(knownPts, knownPos, allTms);
if isnan(allPos(end))
    allPos(end) = allPos(end-1);
end

[~, idx_bin] = max(pxn);
idx_bin = idx_bin.*4;