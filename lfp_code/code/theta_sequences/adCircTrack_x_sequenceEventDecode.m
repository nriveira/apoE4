function adCircTrack_x_sequenceEventDecode(group)
% function fmr1CircTrack_x_sequenceEventDecode(group)
%
% PURPOSE:
%   This code detects sequence events when the rat is running (as in Zheng
%       et al. 2021) and plots the decoded position relative to the rat's
%       current position across the event across normalized event time.
%
% INPUT:
%   group = data struct
%
% OUTPUT:
%   Figures:
%       F1: Decoded position (relative to current) across all sequences
%           from each group.
%       F2: Decoded probability of past, present, and future locations
%           across normalized time in forward events.
%       F3: Decoded probability of past, present, and future locations
%           across normalized time in reverse events.
%       F4: Sequence slopes.
%       F5: Quadrant probability difference from the decoded distribution
%           during sequence event (see Feng et al. 2015 for reference).
%       F6: Correlation between time and decoded position during sequence
%           events AKA weighted correlation.
%       F7: Correlation between spike time and place cell peak firing
%           positions within 50 cm of rat AKA spike train correlation.
%       F8: Sum of the decoded probability distribution (pxn) from current
%           location to reward (see Zheng et al. 2021 for reference).
%
% OPTIONS:
%   See function for all, especially inputs for identifying sequences
%   events. Others of particular note:
%       saveOrNot = 1 to save figs, 0 don't
%       downSampEvents =  whether to down sample cells used in decoding or
%           number of events in each group
%
% MMD
% 7/2021
% Colgin Lab

%% OPTIONS

methToUse = 2; %which method you want to use!
% 1 = from first 3 begins (fourth begin as training set)
% 2 = leave one out method

downSampEvents = 0;

saveOrNot = 1; %for saving main figures ONLY
prepForStats = 0;

bayesStep = 10/1000; %10 ms

plotEachDay = 0;
saveDayPlots = 0; %can make this 1 to save only the day plots and not the main figures

minCellDay = 30; %minimum number of cells in a day to do the decoding

downSampCell = 0; %1 to down sample cells to the same for each day
newCellNum = 40; %number of cells to down sample to

%% INITIALIZE

saveDir = 'E:\resultsFeb2023_AD_WT\thetaSequences';
curDir = pwd;

if downSampEvents == 1 || downSampCell == 1
    rs = RandStream('mt19937ar');
end %get rand stream

spatBinSz = 4;

% normPos = 180; %normalize actual position to be 180
numTimeBins = 100; %normalizing for plot

distBins = 12/spatBinSz; %12 degrees; distance for past and future decoding definitions

degBinCtrs = group(1).rat(1).day(1).binCtrs; %doesn't change across days/rats
radBinCtrs = deg2rad(degBinCtrs);

spatBinSz = 4;
% binMat = 1:length(radBinCtrs);

% degCmConv = 100*pi / 360; %track has 1 m diameter

cols = {'Blue', 'Red'};
seqTypes = {'Forward', 'Reverse'};
groupNames = {'WT', 'FXS'};

radCmConv = (pi*100) / (2*pi); %convert all track measurements from deg to cm before plotting

if plotEachDay == 1 && downSampEvents == 1
    error('Do not down sample events when plotting each day.')
end

for m = methToUse
    fprintf('DECODING USING METHOD %d\n', m)
    cd(saveDir)
    cd(['method' num2str(m)])
    try
        cd(['minCellperDay_' num2str(minCellDay)])
    catch
        mkdir(['minCellperDay_' num2str(minCellDay)])
        cd(['minCellperDay_' num2str(minCellDay)])
    end %make a directory if needed
    if saveDayPlots == 1
        cd('plotEachDay')
    end
    
    %initialze for this method
    
    acSeqPxn = cell(2,2); %group x direction
    probOverEv = cell(2,3,2); %group x past/present/future loc x dirInd
    
    xSpan = cell(2,2); %distances
    xSpanRel = cell(2,2); %x-span - move dist
    tSpan = cell(2,2); %durations
    
    slopeAllSeq = cell(2,2); %n = all sequences for each group! by begin, direction
    
    sumToRew = cell(2,2); %sum of pxn from two bins (10 deg) in front of rat to next reward group x begin x direction
    sumtoLastRew = cell(2,2); %sum of pxn from rat's current pos to last reward
    
    quadProb = cell(2,2); %quandrant sum probability (see Feng et al. 2015)
    weighCorr = cell(2,2); %weighted correlation
    spkTrainCorr = cell(2,2); %spike train correlation
    
    %% PREP FOR DOWN SAMPLE (EVENTS)
    
    if downSampEvents == 1
        minG = zeros(1,2);
        gNum = zeros(2,2);
        keepEvs = cell(2,1);
        for dirInd = 1:2
            for g = 1:2
                for r = 1:length(group(g).rat)
                    for d = 1:length(group(g).rat(r).day)
                        for b = 1:3
                            if ~isempty(group(g).rat(r).day(d).begin(b).seq)
                                gNum(g,dirInd) = gNum(g,dirInd) + size(group(g).rat(r).day(d).begin(b).seq(m).inds{dirInd},1);
                            end %any events
                        end %sleep
                    end %day
                end %rat
            end %group
            
            [minEv, minG(dirInd)] = min(gNum(:,dirInd));
            tmpKeepEvs = datasample(rs, 1:max(gNum(:,dirInd)), minEv, 'Replace', false);
            keepEvs{dirInd} = sort(tmpKeepEvs);
        end %dirInd
        allEvCntr = zeros(1,2); %initialize
    end %if down sampling events
    
    %% GET DATA
    
    for g = 1:2
        fprintf('%s\n', group(g).name)
        for r = 1:length(group(g).rat)
            fprintf('\tRat %d/%d\n', r, length(group(g).rat))
            
            for d = 1:length(group(g).rat(r).day)
                fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))
                
                if isempty(group(g).rat(r).day(d).begin(1).seq)
                    fprintf('\t\t\tNot enough cells to detect sequences\n')
                else
                    fprintf('\t\t\tGetting sequence data...\n')
                    if plotEachDay == 1
                        figtitle = [group(g).rat(r).name '_' group(g).rat(r).day(d).name];
                        figure('Name', figtitle)
                        daySeq = [];
                    end
                    rewLocs = group(g).rat(r).day(d).rewLocs;
                    
                    badU = [];
                    if m == 1
                        unitInfo = group(g).rat(r).day(d).x3BeginUnitInfo;
                    else
                        unitInfo = group(g).rat(r).day(d).xAllBeginUnitInfo;
                    end %which method to get unit info
                    uIDs = zeros(length(unitInfo),2);
                    
                    rateMaps = zeros(length(unitInfo), length(unitInfo(1).smRateMap));
                    uPkPos = zeros(1,length(unitInfo));
                    for u = 1:length(unitInfo)
                        if max(unitInfo(u).smRateMap)>=1 %unit is bad if max firing rate in bin does not exceed 1
                            rateMaps(u,:) = unitInfo(u).smRateMap; %Smoothed ratemap
                            uIDs(u,:) = unitInfo(u).ID;
                            [~, pkInd] = max(unitInfo(u).smRateMap);
                            uPkPos(u) = degBinCtrs(pkInd);
                        else
                            badU = [badU u]; %#ok
                        end %bad fr
                    end %units
                    rateMaps(badU,:) = [];
                    rateMaps(rateMaps==0) = 0.0001; %get rid of zeros because our Bayesian decoder can't handle 'em.
                    uIDs(badU,:) = [];
                    uPkPos(badU) = [];
                    
                    if m == 1
                        begNums = 1:3;
                    else
                        begNums = 1:4;
                    end %which method for begins to cover
                    for b = begNums
                        
                        for dirInd = 1:2
                            fprintf('\t\t\t\tBegin %d/%d\n', b, length(begNums))
                            radPos = group(g).rat(r).day(d).begin(b).radPos;
                            %                         pxn = group(g).rat(r).day(d).begin(b).seq(m).pxn;
                            %
                            %                         seqInds = group(g).rat(r).day(d).begin(b).seq(m).inds;
                            %                         seqTms = group(g).rat(r).day(d).begin(b).seq(m).tms;
                            %                         seqSlopes =  group(g).rat(r).day(d).begin(b).seq(m).slopes;
                            
                            
                            if downSampEvents == 1 && g ~= minG(dirInd) %find out which events for this direction from this begin to keep
                                keyboard
                                begNum = size(seqInds{dirInd},1); %number of events in this begin
                                begEvs = 1:begNum;
                                keepInds = keepEvs{dirInd} > allEvCntr(dirInd) & keepEvs{dirInd} <= (allEvCntr(dirInd) + begNum);
                                
                                begEvInds = keepEvs{dirInd}(keepInds) - allEvCntr(dirInd);
                                begEvs = begEvs(begEvInds);
                                
                                allEvCntr(dirInd) = allEvCntr(dirInd) + begNum;
                            end %find out fr
                            
                            
                            seq = group(g).rat(r).day(d).begin(b).method(m).seqDir(dirInd).seq;
                            
                            begEvCntr = 0;
                            for i = 1:length(seq)
                                begEvCntr = begEvCntr + 1;
                                
                                
                                if downSampEvents == 0 || g == minG(dirInd) || g ~= minG(dirInd) && ismember(begEvCntr, begEvs)
                                    
                                    %                                 sqOnInd = seqInds{dirInd}(i,1);
                                    %                                 sqOffInd = seqInds{dirInd}(i,2);
                                    
                                    seqPxn = seq(i).pxn;
                                    %                                 seqPxn(isnan(seqPxn)) =
                                    %                                 1/size(rateMaps,2); %change nans to just
                                    %                                 chance - should be none
                                    pxnTAxis = bayesStep:bayesStep:(diff(seq(i).inds)+1)*bayesStep;
                                    
                                    %                                 = seqTms{dirInd}(i,1); %get actual time
                                    %                                 startInd = match(sqOnTm, radPos(:,1));
                                    %                                 if isempty(startInd)
                                    %                                     startInd = 1;
                                    %                                 end
                                    %                               sqOnTm = seq(i).tms(1);
                                    %                                 sqOffTm = seq(i).tms(2);
                                    
                                    if diff(seq(i).tms) > 0.15
                                        continue
                                    end
                                    
                                    %                                 endInd = match(sqOffTm, radPos(:,1));
                                    %                                 actPos = mean(radPos(startInd:endInd,2)); %get position of animal during sequence event
                                    actInd = round(seq(i).actPos/deg2rad(spatBinSz));
                                    
                                    tmpxSpan = abs(circ_dist(seq(i).com(end), seq(i).com(1)));
                                    tmptSpan = diff(seq(i).tms);
                                    
                                    startInd = match(seq(i).tms(1), radPos(:,1));
                                    endInd = match(seq(i).tms(2), radPos(:,1));
                                    
                                    moveDist = abs(circ_dist(deg2rad(radPos(endInd,2)), deg2rad(radPos(startInd,2))));
                                    tmpxSpanRel = tmpxSpan - moveDist;
                                    
                                    xSpan{g,dirInd} = [xSpan{g,dirInd} tmpxSpan];
                                    xSpanRel{g,dirInd} = [xSpanRel{g,dirInd} tmpxSpanRel];
                                    tSpan{g,dirInd} = [tSpan{g,dirInd} tmptSpan];
                                    
                                    slopeAllSeq{g,dirInd} = [slopeAllSeq{g,dirInd}; abs(seq(i).slope)*radCmConv];
                                    
                                    distToRew = rad2deg(circ_dist(deg2rad(rewLocs), seq(i).actPos));
                                    if length(rewLocs) == 1 && distToRew < 0
                                        tmpSumToRew = []; %only do it within 180 deg for one reward zone, to be consistent
                                    else
                                        tmpSumToRew = seqSumPxn_toRew(seqPxn, rewLocs, rad2deg(seq(i).actPos));
                                        sumToRew{g,dirInd} = [sumToRew{g,dirInd} tmpSumToRew];
                                    end %reward zone check
                                    
                                    shPxn = shift_pxn(seqPxn, rad2deg(seq(i).actPos), spatBinSz);
                                    
                                    tmpPxn = zeros(length(radBinCtrs),numTimeBins); %rescale for normalize time in sequence event - since events can be dif times
                                    tmpBinTm = rescale(pxnTAxis, 0, 1);
                                    
                                    newBinStart = 1;
                                    for bn = 2:length(tmpBinTm)
                                        newBinEnd = round(tmpBinTm(bn),2);
                                        newBinEnd = round(newBinEnd/(1/numTimeBins),0);
                                        
                                        tmpPxn(:,newBinStart:newBinEnd) = repmat(shPxn(:,bn-1), [1, newBinEnd-newBinStart+1]); %tile to fit norm time distribution
                                        newBinStart = newBinEnd;
                                    end %bin
                                    
                                    acSeqPxn{g,dirInd} = cat(3, acSeqPxn{g,dirInd}, tmpPxn);
                                    
                                    sumPast = sum(tmpPxn(size(tmpPxn,1)/2-distBins-1:size(tmpPxn,1)/2-distBins+1,:));
                                    sumPres = sum(tmpPxn(size(tmpPxn,1)/2-(distBins-1)/2:size(tmpPxn,1)/2+(distBins-1)/2,:));
                                    sumFut = sum(tmpPxn(size(tmpPxn,1)/2+distBins-1:size(tmpPxn,1)/2+distBins+1,:));
                                    
                                    %                                     sumPast = sum(shPxn(size(shPxn,1)/2-distBins-1:size(shPxn,1)/2-distBins+1,:));
                                    %                                     sumPres = sum(shPxn(size(shPxn,1)/2-(distBins-1)/2:size(shPxn,1)/2+(distBins-1)/2,:));
                                    %                                     sumFut = sum(shPxn(size(shPxn,1)/2+distBins-1:size(shPxn,1)/2+distBins+1,:));
                                    
                                    probOverEv{g,1,dirInd} = [probOverEv{g,1,dirInd}; sumPast];
                                    probOverEv{g,2,dirInd} = [probOverEv{g,2,dirInd}; sumPres];
                                    probOverEv{g,3,dirInd} = [probOverEv{g,3,dirInd}; sumFut];
                                    
                                    if plotEachDay == 1 && dirInd == 1
                                        daySeq = cat(3, daySeq, tmpPxn);
                                    end
                                    probDiff = seqProbabilityDiff(shPxn, spatBinSz);
                                    quadProb{g,dirInd} = [quadProb{g,dirInd} probDiff];
                                    %For visualizing 4 quadrants (see Feng et al. 2015)
                                    %imagesc(1:size(shiftPxn,2),-180:180,shiftPxn); axis xy
                                    %line([1 size(shiftPxn,2)], [0 0], 'LineStyle', '--', 'Color', 'white')
                                    %line([ceil(size(shiftPxn,2)/2) ceil(size(shiftPxn,2)/2)], [-180 180], 'LineStyle', '--', 'Color', 'white')
                                    %ylim([-ceil(degCmConv*50) ceil(degCmConv*50)])
                                    
                                    corrTP = seqWeighCorr(shPxn, spatBinSz);
                                    weighCorr{g,dirInd} = [weighCorr{g,dirInd} corrTP];
                                    
                                    %find place cells with pk pos within 50 cm of act
                                    %                                 shiftUPkPos = wrapTo360(uPkPos + shiftDeg);
                                    %                                 spkTmsForCorr = [];
                                    %                                 pkPosForCorr = [];
                                    seqSpkTms = cell(1,length(uIDs));
                                    uCntr = 0;
                                    for u = 1:length(group(g).rat(r).day(d).x3BeginUnitInfo)
                                        uID = group(g).rat(r).day(d).x3BeginUnitInfo(u).ID;
                                        if ismember(uID, uIDs, 'row')
                                            uCntr = uCntr + 1;
                                            seqSpks = group(g).rat(r).day(d).begin(b).unit(u).spkTms(group(g).rat(r).day(d).begin(b).unit(u).spkTms >= seq(i).tms(1) & group(g).rat(r).day(d).begin(b).unit(u).spkTms <= seq(i).tms(2));
                                            
                                            seqSpkTms{uCntr} = seqSpks;
                                        end %use unit
                                    end %unit
                                    corrSP = seqSpikeTrainCorr(seqSpkTms, uPkPos, rad2deg(seq(i).actPos));
                                    
                                    if ~isnan(corrSP)
                                        spkTrainCorr{g,dirInd} = [spkTrainCorr{g,dirInd} corrSP];
                                    end
                                    
                                end %if event should be included (when down sampling)
                            end %sequences
                        end %dirInd
                    end %begin
                    
                    if plotEachDay == 1
                        dayMean = mean(daySeq,3);
                        imagesc(0:1, -180:spatBinSz:180, dayMean)
                        axis xy
                        axis square
                        
                        [~, maxInds] = max(dayMean);
                        
                        pFit = polyfit(1/100:1/100:1, degBinCtrs(maxInds), 1);
                        hold on;
                        yFit = [pFit(2) pFit(1)*1+pFit(2)] - 180; %since 0 is now current location
                        plot([0 1], yFit, 'w')
                        
                        ylim([-55 55])
                        %line([0 1], [0 0], 'Color', rgb('white'), 'LineWidth', 2, 'LineStyle', '--')
                        yticks(-90:45:90)
                        ylabel('Position relative to actual (degrees)')
                        
                        xlabel('Normalized time in sequence event')
                        xticks([0 1])
                        
                        colormap(jet)
                        cbr = colorbar;
                        ylabel(cbr, 'Probability')
                        title([group(g).name ' Rat ' num2str(r) ' Day ' num2str(d) ' (n = ' num2str(size(daySeq,3)) ' seq)'])
                        
                        if saveDayPlots == 1
                            saveas(gcf, figtitle, 'epsc')
                            saveas(gcf, figtitle, 'png')
                            saveas(gcf, figtitle, 'fig')
                        end
                    end %plot each day
                    
                end %begin
            end %day
        end %rat
    end %group
    
    if saveOrNot == 1 && saveDayPlots == 1
        cd ../
    elseif saveOrNot == 0 && saveDayPlots == 1
        cd(curDir)
    end %get to our correct directory
    keyboard
    %% FIG 1 - PXN ACR SEQ
    
    figtitle = 'DecodeAcrossSequenceEvents_byEvent';
    if downSampCell == 1
        figtitle = [figtitle '_downSample_' num2str(newCellNum) 'cellsPerDay'];
    end
    if downSampEvents == 1
        figtitle = [figtitle '_downSampled_events'];
    end %down sample events
    
    tmpMax = 0;
    for g = 1:2
        for dirInd = 1:2
            if max(max(mean(acSeqPxn{g,dirInd},3))) > tmpMax
                tmpMax =  max(max(mean(acSeqPxn{g,dirInd},3)));
            end
        end %dirInd
    end %group
    
    figure('Name', figtitle, 'Position', [336 160 1378 753])
    spMap = [1 3; 2 4];
    for g = 1:2
        for dirInd = 1:2
            subplot(2,2,spMap(g,dirInd))
            
            imagesc(0:1, -180:spatBinSz:180, mean(acSeqPxn{g,dirInd},3))
            axis xy
            axis square
            
            %             [~, maxInds] = max(mean(acSeqPxn{g,dirInd},3));
            %
            %             pFit = polyfit(1/100:1/100:1, degBinCtrs(maxInds), 1);
            %             hold on;
            %             yFit = [pFit(2) pFit(1)*1+pFit(2)] - 180; %since 0 is now current location
            %             plot([0 1], yFit, 'w', 'LineWidth', 2)
            
            ylim([-65 65])
            yticks(-90:45:90)
            %line([0 1], [0 0], 'Color', rgb('white'), 'LineWidth', 2, 'LineStyle', '--')
            
            if g == 1
                ylabel({[seqTypes{dirInd} ' Sequences'], 'Position relative to actual (degrees)'})
            else
                ylabel('Position relative to actual (degrees)')
            end
            xlabel('Normalized time in sequence event')
            xticks([0 1])
            
            colormap(jet)
            cbr = colorbar;
            ylabel(cbr, 'Probability')
            caxis([0 tmpMax])
            cbr.Ticks = 0:0.02:tmpMax;
            
            title([group(g).name ' n = ' num2str(size(acSeqPxn{g,dirInd},3)) ' sequences'])
            
        end %slopeInd
    end %group
    
    if saveOrNot == 1
        saveas(gcf, figtitle, 'epsc')
        saveas(gcf, figtitle, 'png')
        saveas(gcf, figtitle, 'fig')
    end %save option
    
    %% FIG 2 - RECONSTRUCT DECODE ACROSS EVENT
    
    figtitle = 'DecodeAcrossSequenceEvents_byEvent_reconstructed';
    if downSampEvents == 1
        figtitle = [figtitle '_downSampled_events'];
    end %down sample events
    
    rotBins = deg2rad(-180+spatBinSz/2:spatBinSz:180); %rotate bins
    
    figure('Name', figtitle, 'Position', [336 160 1378 753])
    spMap = [1 3; 2 4];
    for g = 1:2
        for dirInd = 1:2
            subplot(2,2,spMap(g,dirInd))
            
            tmpMean = mean(acSeqPxn{g,dirInd},3);
            com = zeros(1, size(tmpMean,2));
            for t = 1:size(tmpMean,2)
                com(t) = rad2deg(circ_mean(rotBins', tmpMean(:,t)));
            end %time bin
            
            plot(linspace(0,1,length(com)), com, 'k.')
            axis square
            
            beta = CircularRegression(linspace(0,1,length(com)), deg2rad(com));
            cmSlope = beta(1) * radCmConv;
            calphase = beta(1)*linspace(0,1,length(com)) + beta(2);
            calphase = rad2deg(calphase);
            
            hold on;
            plot(linspace(0,1,length(com)), calphase, 'b')
            
            ylim([-25 25])
            
            if g == 1
                ylabel({[seqTypes{dirInd} ' Sequences'], 'Position relative to actual (degrees)'})
            else
                ylabel('Position relative to actual (degrees)')
            end
            xlabel('Normalized time in sequence event')
            xticks([0 1])
            
            title({[group(g).name ' n = ' num2str(size(acSeqPxn{g,dirInd},3)) ' sequences'], ['slope = ' num2str(round(cmSlope)) ' cm/s']})
            
        end %slopeInd
    end %group
    
    if saveOrNot == 1
        saveas(gcf, figtitle, 'epsc')
        saveas(gcf, figtitle, 'png')
        saveas(gcf, figtitle, 'fig')
    end %save option
    
%     %% FIG - X-SPAN
%     
%     figtitle = 'Seq_xSpan_distance';
%     figure('Name', figtitle, 'Position', [398 318 1028 648])
%     
%     meanData = cellfun(@mean, xSpan);
%     semData = cellfun(@semfunct, xSpan);
%     
%     for dirInd = 1:2
%         subplot(2,2,dirInd)
%         
%         bgraph = bar(meanData(:,dirInd), 'FaceColor', 'Flat');
%         hold on;
%         errorbar(1:2, meanData(:,dirInd), semData(:,dirInd), semData(:,dirInd), 'LineStyle', 'None', 'Color', 'Black')
%         for g = 1:2
%             bgraph.CData(g,:) = rgb(cols{g});
%         end %group
%         
%         ylabel('x-span (rad)')
%         xticklabels(groupNames)
%         title([seqTypes{dirInd} ' sequences'])
%     end %dirInd
%     
%      meanData = cellfun(@mean, xSpanRel);
%     semData = cellfun(@semfunct, xSpanRel);
%     
%     for dirInd = 1:2
%         subplot(2,2,dirInd+2)
%         
%         bgraph = bar(meanData(:,dirInd), 'FaceColor', 'Flat');
%         hold on;
%         errorbar(1:2, meanData(:,dirInd), semData(:,dirInd), semData(:,dirInd), 'LineStyle', 'None', 'Color', 'Black')
%         for g = 1:2
%             bgraph.CData(g,:) = rgb(cols{g});
%         end %group
%         
%         ylabel('Relative x-span (rad)')
%         xticklabels(groupNames)
%         title([seqTypes{dirInd} ' sequences'])
%     end %dirInd
%     
%     if saveOrNot == 1
%             saveas(gcf, figtitle, 'epsc')
%             saveas(gcf, figtitle, 'png')
%             saveas(gcf, figtitle, 'fig')
%         end %save option
%         
%     %% FIG - T SPAN
%     
%      figtitle = 'Seq_tSpan_duration';
%     figure('Name', figtitle, 'Position', [397 640 1028 327])
%     
%     meanData = cellfun(@mean, tSpan);
%     semData = cellfun(@semfunct, tSpan);
%     
%     for dirInd = 1:2
%         subplot(1,2,dirInd)
%         
%         bgraph = bar(meanData(:,dirInd), 'FaceColor', 'Flat');
%         hold on;
%         errorbar(1:2, meanData(:,dirInd), semData(:,dirInd), semData(:,dirInd), 'LineStyle', 'None', 'Color', 'Black')
%         for g = 1:2
%             bgraph.CData(g,:) = rgb(cols{g});
%         end %group
%         
%         ylabel('t-span (s)')
%         xticklabels(groupNames)
%         title([seqTypes{dirInd} ' sequences'])
%     end %dirInd
%     
%        if saveOrNot == 1
%             saveas(gcf, figtitle, 'epsc')
%             saveas(gcf, figtitle, 'png')
%             saveas(gcf, figtitle, 'fig')
%         end %save option
%     
%     %% FIG 3 - PROBILITY
%     
%     for dirInd = 1:2
%         figtitle = ['Decoded_pastPresFut_acrossEvents_' seqTypes{dirInd}];
%         figure('Name', figtitle, 'Position', [413 637 1007 326])
%         
%         tmpCols = {'Orange', 'Green', 'Purple'};
%         decNames = {'Past' 'Present', 'Future'};
%         nboot = 5000;
%         for g = 1:2
%             subplot(1,2,g)
%             lh = zeros(1,3);
%             for t = 1:3
%                 tmpData = probOverEv{g,t,dirInd};
%                 meanData = mean(tmpData);
%                 CI = bootci(nboot, {@mean, tmpData}, 'type', 'per');
%                 
%                 lh(t)=  plot_filled_ci(linspace(0,1,100), meanData, CI, rgb(tmpCols{t}));
%             end %decoded types
%             
%             ylabel('Decoded probability')
%             xlabel('Normalized time in sequence event')
%             title([group(g).name ' (n = ' num2str(size(probOverEv{g,1,dirInd},1)) ' seq)'])
%             axis square;
%         end %group
%         same_axes;
%         leg = legend(lh, decNames, 'Location', 'Northeastoutside');
%         leg.Position = [0.93 0.75 0.02 0.1];
%         
%         if saveOrNot == 1
%             saveas(gcf, figtitle, 'epsc')
%             saveas(gcf, figtitle, 'png')
%             saveas(gcf, figtitle, 'fig')
%         end %save option
%     end %dir ind
%     
%     %% FIG 3 - SLOPES
%     
%     %         figtitle = 'SequenceSlopes';
%     %         if downSampEvents == 1
%     %             figtitle = [figtitle '_downSampEvents'];
%     %         end
%     %
%     %         figure('Name', figtitle, 'Position', [478 644 793 321])
%     %
%     %         yLab = 'Slope (|cm/s|)';
%     %         xLab = 'Event type';
%     %         group1Names = {'WT', 'FXS'};
%     %          group2Names = {'Forward', 'Reverse'};
%     %         bar_and_dotplot_grouped(slopeAllSeq, yLab, xLab, group1Names, group2Names, cols)
%     %
%     
%     %     meanSlope = rad2deg(cellfun(@mean, slopeAllSeq));
%     meanSlope = cellfun(@mean, slopeAllSeq);
%     errorSlope = cellfun(@semfunct, slopeAllSeq);
%     
%     figtitle = 'SequenceSlopes';
%     if downSampEvents == 1
%         figtitle = [figtitle '_downSampEvents'];
%     end
%     
%     figure('Name', figtitle, 'Position', [478 644 793 321])
%     
%     for dirInd = 1:2
%         subplot(1,2,dirInd)
%         
%         bgraph = bar(meanSlope(:,dirInd), 'FaceColor', 'Flat');
%         hold on;
%         errorbar(1:2, meanSlope(:,dirInd), errorSlope(:,dirInd), errorSlope(:,dirInd), 'LineStyle', 'None', 'Color', 'Black')
%         
%         xLabs = cell(1,2);
%         for g = 1:2
%             bgraph.CData(g,:) = rgb(cols{g});
%             
%             xLabs{g} = [group(g).name ' n = ' num2str(length(slopeAllSeq{g,dirInd}))];
%         end
%         
%         ylabel('Slope (|cm/s|)')
%         %         ylabel('Slope (|rad/s|)')
%         xticklabels(xLabs)
%         title([seqTypes{dirInd} ' Sequences'])
%     end %dirInd
%     
%     if saveOrNot == 1
%         saveas(gcf, figtitle, 'epsc')
%         saveas(gcf, figtitle, 'png')
%         saveas(gcf, figtitle, 'fig')
%     end %save option
%     
%     %% FIG 4 - SLOPES DISTRIBUTIONS - FORWARD ONLY
%     
%     figtitle = 'SequenceSlopes_distributions';
%     figure('Name', figtitle, 'Position', [478 644 793 321])
%     dirInd = 1;
%     
%     minSlope = floor(min(vertcat(slopeAllSeq{:,dirInd})));
%     maxSlope = ceil(max(vertcat(slopeAllSeq{:,dirInd})));
%     
%     subplot(1,2,1)
%     leg = cell(2,1);
%     lh = nan(2,1);
%     for g = 1:2
%         [cumProp, xScale] = calc_cum_prop_ip_range(slopeAllSeq{g,dirInd}, [minSlope maxSlope]);
%         lh(g) = plot(xScale, cumProp, 'Color', rgb(cols{g}));
%         leg{g} = [groupNames{g} ' ( n = ' num2str(length(slopeAllSeq{g,dirInd})) ' seq)'];
%         hold on;
%         
%         line([meanSlope(g,dirInd) meanSlope(g,dirInd)], [0 1], 'Color', rgb(cols{g}), 'LineStyle', '--')
%     end %group
%     
%     ylabel('Proportion')
%     xlabel('Slope (cm/s')
%     legend(lh, leg, 'Location', 'southeast')
%     
%     subplot(1,2,2)
%     sigma = 50;
%     for g = 1:2
%         [weighDis, xVals] = WeightedProportion(slopeAllSeq{g,dirInd}, minSlope, maxSlope, sigma); %from Ernie's code for Zheng et al. 2021
%         plot(xVals, weighDis, 'Color', rgb(cols{g}))
%         hold on;
%     end %group
%     
%     ax = gca;
%     yMax = ax.YLim(2);
%     for g = 1:2
%         line([meanSlope(g,dirInd) meanSlope(g,dirInd)], [0 yMax], 'Color', rgb(cols{g}), 'LineStyle', '--')
%     end %group
%     
%     ylabel('Proportion')
%     xlabel('Slope (cm/s)')
%     xlim([0 maxSlope])
%     legend(leg, 'Location', 'northeast')
%     
%     if saveOrNot == 1
%         saveas(gcf, figtitle, 'epsc')
%         saveas(gcf, figtitle, 'png')
%         saveas(gcf, figtitle, 'fig')
%     end %save option
%     
%     %% FIG 4 - PROBABILITY DIFFERENCE
%     
%     figtitle = 'ProbabilityDifferences';
%     if downSampEvents == 1
%         figtitle = [figtitle '_downSampEvents'];
%     end
%     
%     figure('Name', figtitle, 'Position', [478 644 793 321])
%     
%     meanProb = cellfun(@mean, quadProb);
%     semProb = cellfun(@semfunct, quadProb);
%     
%     for dirInd = 1:2
%         subplot(1,2,dirInd)
%         
%         bgraph = bar(meanProb(:,dirInd), 'FaceColor', 'Flat');
%         
%         hold on;
%         er = errorbar(1:2, meanProb(:,dirInd), semProb(:,dirInd), semProb(:,dirInd), 'Color', 'Black', 'LineStyle', 'None');
%         
%         xLabs = cell(1,2);
%         for g = 1:2
%             bgraph.CData(g,:) = rgb(cols{g});
%             xLabs{g} = [group(g).name ' n = ' num2str(length(quadProb{g,dirInd}))];
%         end
%         
%         xticklabels(xLabs)
%         ylabel('(II+IV-I-III)/sum quadrant probability')
%         title([seqTypes{dirInd} ' Sequences'])
%         
%     end %dirInd
%     
%     if saveOrNot == 1
%         saveas(gcf, figtitle, 'epsc')
%         saveas(gcf, figtitle, 'png')
%         saveas(gcf, figtitle, 'fig')
%     end %save option
%     
%     %% FIG 5 - WEIGHTED CORRELATION
%     
%     figtitle = 'WeightedCorrelation';
%     if downSampEvents == 1
%         figtitle = [figtitle '_downSampEvents'];
%     end
%     
%     figure('Name', figtitle, 'Position', [478 644 793 321])
%     
%     meanCorr = cellfun(@mean, weighCorr);
%     semCorr = cellfun(@semfunct, weighCorr);
%     
%     for dirInd = 1:2
%         subplot(1,2,dirInd)
%         
%         bgraph = bar(meanCorr(:,dirInd), 'FaceColor', 'Flat');
%         
%         hold on;
%         er = errorbar(1:2, meanCorr(:,dirInd), semCorr(:,dirInd), semCorr(:,dirInd), 'Color', 'Black', 'LineStyle', 'None');
%         
%         xLabs = cell(1,2);
%         for g = 1:2
%             bgraph.CData(g,:) = rgb(cols{g});
%             xLabs{g} = [group(g).name ' n = ' num2str(length(weighCorr{g,dirInd}))];
%         end
%         
%         xticklabels(xLabs)
%         ylabel('Weighted correlation')
%         title([seqTypes{dirInd} ' Sequences'])
%     end %dirInd
%     
%     if saveOrNot == 1
%         saveas(gcf, figtitle, 'epsc')
%         saveas(gcf, figtitle, 'png')
%         saveas(gcf, figtitle, 'fig')
%     end %save option
%     
%     %% FIG 6 - SPIKE TRAIN CORRELATION
%     
%     figtitle = 'SpikeTrainCorrelation';
%     if downSampEvents == 1
%         figtitle = [figtitle '_downSampEvents'];
%     end
%     
%     figure('Name', figtitle, 'Position', [478 644 793 321])
%     
%     meanCorr = cellfun(@nanmean, spkTrainCorr);
%     semCorr = cellfun(@nansemfunct, spkTrainCorr);
%     
%     for dirInd = 1:2
%         subplot(1,2,dirInd)
%         
%         bgraph = bar(meanCorr(:,dirInd), 'FaceColor', 'Flat');
%         
%         hold on;
%         er = errorbar(1:2, meanCorr(:,dirInd), semCorr(:,dirInd), semCorr(:,dirInd), 'Color', 'Black', 'LineStyle', 'None');
%         
%         xLabs = cell(1,2);
%         for g = 1:2
%             bgraph.CData(g,:) = rgb(cols{g});
%             xLabs{g} = [group(g).name ' n = ' num2str(length(spkTrainCorr{g,dirInd}))];
%         end
%         
%         xticklabels(xLabs)
%         ylabel('Spike train correlation')
%         title([seqTypes{dirInd} ' Sequences'])
%         
%     end %dirInd
%     
%     if saveOrNot == 1
%         saveas(gcf, figtitle, 'epsc')
%         saveas(gcf, figtitle, 'png')
%         saveas(gcf, figtitle, 'fig')
%     end %save option
%     
%     %% FIG 7 - SUM(PXN) TO REWARD
%     
%     figtitle = 'SumPxntoRew';
%     if downSampEvents == 1
%         figtitle = [figtitle '_downSampEvents'];
%     end
%     
%     figure('Name', figtitle, 'Position', [478 644 793 321])
%     
%     meanSumToRew = cellfun(@mean, sumToRew);
%     semSumToRew = cellfun(@semfunct, sumToRew);
%     
%     for dirInd = 1:2
%         subplot(1,2,dirInd)
%         
%         bgraph = bar(meanSumToRew(:,dirInd), 'FaceColor', 'Flat');
%         
%         hold on;
%         er = errorbar(1:2, meanSumToRew(:,dirInd), semSumToRew(:,dirInd), semSumToRew(:,dirInd), 'Color', 'Black', 'LineStyle', 'None');
%         
%         xLabs = cell(1,2);
%         for g = 1:2
%             bgraph.CData(g,:) = rgb(cols{g});
%             xLabs{g} = [group(g).name ' n = ' num2str(length(weighCorr{g,dirInd}))];
%         end
%         
%         xticklabels(xLabs)
%         ylabel('Sum(Pxn) to reward')
%         title([seqTypes{dirInd} ' Sequences'])
%     end %dirInd
%     
%     if saveOrNot == 1
%         saveas(gcf, figtitle, 'epsc')
%         saveas(gcf, figtitle, 'png')
%         saveas(gcf, figtitle, 'fig')
%     end %save option
%     
%     
%     %% STATS?
%     
%     if prepForStats == 1
%         keyboard
%         statSlopes = cell(2,1);
%         for dirInd = 1:2
%             for g = 1:2
%                 
%                 for i = 1:length(slopeAllSeq{g,dirInd})
%                     statSlopes{dirInd} = [statSlopes{dirInd}; g  slopeAllSeq{g,dirInd}(i)];
%                 end %sequences
%                 
%             end %group
%         end %dirInd
%         keyboard
%     end %prep for stats
    
end %method to use

cd(curDir)

end %function