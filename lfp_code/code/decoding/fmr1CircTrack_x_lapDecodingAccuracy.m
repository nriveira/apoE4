function fmr1CircTrack_x_lapDecodingAccuracy(group)
% function fmr1CircTrack_x_lapDecodingAccuracy(group)
%
% *See NOTE below
%
% PURPOSE:
%  Function generates a cumulative distribution plot of decoding error with a separate line for each
%  day/rat, with separate subplots for each group (WT vs KO).
%  - Function optionally also plots most probable location based on decoding on top of PPM heat mats
%    with rat's actual position overlaid on that (plotActVsDecPosnEachLap)
%  - Function also optionally plots error averaged across laps for each 'begin'
%
% INPUT:
%  group = project uber data struct
%  - Internal options at top of function allow for plotting/saving of optional graphs described
%    within PURPOSE
%
% OUTPUT:
%  See 'PURPOSE' above.
%
%  NOTE!!!
%   This function can not plot a 'correct' confusion matrix because some units have ratemaps showing peaks
%   around reward 1 (0/360*), but laps are defined as when the rat exits/enters +/-15* from 0/360. Therefore
%   the rat's actual position is never within this 345:+15* space, but the decoder will say the rat *is* there
%   and that leads to significant issues that ultimately make the confusion matrices have very high probabilities
%   near 0/360, which is artificial. SO! use 'fmr1CircTrack_x_confusionMatrix' instead. This option is kept here
%   despite these issues for the sake of comparison.
%
% JB Trimper
% 01/2021
% Colgin Lab



% Plot PPM, decoded position, and actual position for each lap (1 subplot per lap for each 'begin')
%   + save options
plotActVsDecPosnEachLap = 0;
saveLapPlots = 0;
saveDir1 = 'C:\Users\John\Desktop\LAB_STUFF\PROJECTS\FMR1_CircTrack\RESULTS\decodedLapFigs';

% Plot error averaged across laps for each 'begin'
%   + save options
plotErrByBegin = 0;
saveBeginErrPlots = 0;
saveDir2 = 'C:\Users\John\Desktop\LAB_STUFF\PROJECTS\FMR1_CircTrack\RESULTS\decodeAccuracyByBegin';

%Plot confusion matrix for each day
%   + save options
plotConfMatsByDay = 0; %Set to 1 to plot confusion matrices for each day, with cell counts in the title
saveDayPlots = 0; %Set to 1 to save these plots automatically in...
saveDir3 = 'C:\Users\John\Desktop\LAB_STUFF\PROJECTS\FMR1_CircTrack\RESULTS\confusionMatricesByDay_byLapMethod';

%Options for saving the across session confusion matrix plots
saveAllSesConfMats = 0; %Set to 1 to save confusion matrix where data is concatenated across all days
saveDir4 = 'C:\Users\John\Desktop\LAB_STUFF\PROJECTS\FMR1_CircTrack\RESULTS\confusionMatrices_acrossAllSes';

% Options for saving the cumulative distribution plot
saveCumErrPlots = 0;
saveDir5 = 'C:\Users\John\Desktop\LAB_STUFF\PROJECTS\FMR1_CircTrack\RESULTS\decodeCumulativeErrAcrossSessions';

% Decoding parameters
bayesWin = 0.50; %s
bayesStep = 0.1; %ms
sampRate = 2000; %Hz
runThresh = 5; %cm/s

radBinCtrs = group(2).rat(1).day(1).binCtrs; %doesn't change across days/rats
newRewLoc = 0; %used to get all plots to line up as though rewards are at 0 & 180
[~,newRewInd] = min(abs(circ_dist(deg2rad(radBinCtrs), deg2rad(newRewLoc))-0));

dayDecErrFig = figure('name', 'Day Decoding Error', 'Position', [25 207.4 1220.8 506.4]);
confMatFig = figure('Name', 'Confusion Matrix', 'Position', [320 248 1124 511]);
for g = 1:2
    fprintf('Group %d\n', g);
    
    figure(dayDecErrFig);
    subplot(1,2,g);
    xlim([0 180]);
    xlabel('Decoding Error (degrees)');
    if g == 1
        ylabel('Cumulative Proportion');
    end
    title(group(g).name);
    fix_font;
    axis square;
    hold on;
    
    allConfMats = [];
    
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        
        
        dayNums = 1:length(group(g).rat(r).day);
        if g == 2
            dayNums = 2; %So we're only considering the KO rat's day w/ most cells
        end
        for d = dayNums
            %         for d = 1:length(group(g).rat(r).day);
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            
            
            confMat = zeros(length(radBinCtrs),length(radBinCtrs)); %Confusion Matrix sums
            posCount = zeros(length(radBinCtrs),1); %For the denominator to calc confMatrix probability
            dayDecodeErr = [];
            
            % Get ratemaps, removing ratemaps if units did not have peak fr >= 1 Hz
            rateMaps = zeros(length(group(g).rat(r).day(d).xAllBeginUnitInfo), length(group(g).rat(r).day(d).xAllBeginUnitInfo(1).smRateMap));
            badU = [];
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                if max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).rateMap)>=1
                    rateMaps(u,:) = group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap; %Smoothed ratemap
                else
                    badU = [badU u]; %#ok
                end
            end
            rateMaps(badU,:) = [];
            rateMaps(rateMaps==0) = 0.0001; %get rid of zeros because our Bayesian decoder can't handle 'em.
            
            if size(rateMaps,1) > 1 %If there are at least 2 units...
                
                % Shift the rate maps so, across all data, reward locations are always plotted at 0 & 180
                %   So plots here will always look like rat goes from 0 to 360
                rewLoc = group(g).rat(r).day(d).rewLocs(1);
                [~,rewInd] = min(abs(circ_dist(deg2rad(radBinCtrs), deg2rad(rewLoc))-0));
                shiftVal = newRewInd - rewInd;
                rateMaps = circshift(rateMaps, [0 shiftVal]);
                
                if plotErrByBegin == 1
                    %             errFig = figure('Position', [ 65   734   972   244]);
                    errFig = figure('name', 'Error Figure', 'Position', [66.6   311.4   1150.4   316]);
                    yMaxes = zeros(1,4);
                end
                
                
                for b = 1:4
                    fprintf('\t\t\tBegin %d\n', b)
                    
                    if plotActVsDecPosnEachLap == 1
                        lapFig = figure('name', 'Lap Figure', 'Position', [ 1    247   1536  757], 'Name', ['Begin ' num2str(b)]);
                    end
                    
                    
                    %Get rat's speed
                    radPos = group(g).rat(r).day(d).begin(b).radPos;
                    instRs = get_runspeed(group(g).rat(r).day(d).begin(b).coords);
                    smRs = smooth_runspeed(instRs);
                    
                    avgErr = nan(1,size(group(g).rat(r).day(d).begin(b).lapInds,1));
                    semErr = nan(1,size(group(g).rat(r).day(d).begin(b).lapInds,1));
                    for lp = 1:size(group(g).rat(r).day(d).begin(b).lapInds,1)
                        
                        startLapTm = group(g).rat(r).day(d).begin(b).lapTms(lp,1);
                        endLapTm = group(g).rat(r).day(d).begin(b).lapTms(lp,end);
                        
                        spkRstr = zeros(length(group(g).rat(r).day(d).begin(b).unit), round((endLapTm-startLapTm)*sampRate));
                        
                        for u = 1:length(group(g).rat(r).day(d).begin(b).unit)
                            allSpkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                            lapSpks = allSpkTms(allSpkTms>=startLapTm & allSpkTms<=endLapTm);
                            timePassed = lapSpks - startLapTm;
                            spkInds = round(timePassed * sampRate);
                            spkInds(spkInds==0) = []; %in case a spike occurred before the first video time-stamp
                            spkRstr(u,spkInds) = 1;
                        end
                        spkRstr(badU,:) = []; %Get rid of units w/ FR that is too low
                        
                        
                        % Bayesian Decoder
                        ppm = BayesianDecoder(spkRstr, rateMaps, bayesWin, bayesStep, sampRate);
                        ppm(isnan(ppm)) = 1 / size(rateMaps,2);
                        
                        % Figure out the number of time windows, their start and middle times
                        [nWin, winStartInds] = find_num_windows(size(spkRstr,2), bayesWin*sampRate, bayesStep*sampRate);
                        if winStartInds(end)+bayesWin*sampRate < size(spkRstr,2)
                            nWin = nWin + 1;
                            winStartInds(end+1) = winStartInds(end)+bayesStep*sampRate; %#ok
                        end
                        winStartTms = startLapTm + winStartInds/sampRate;
                        winEndTms = winStartTms + bayesWin;
                        
                        winMidTms = winStartTms + bayesWin / 2;
                        if winEndTms(end)> endLapTm
                            winEndTms(end) = endLapTm;
                            winMidTms(end) = endLapTm; %because the last time window can extend past the actual endApp time
                        end
                        
                        
                        % Get the rat's actual position
                        %  Don't give a position for the time window if the rat wasn't moving > runThresh
                        actPosn = nan(1,length(winStartTms));
                        for i = 1:length(winStartTms)
                            winSpd = mean(smRs(smRs(:,1)>=winStartTms(i) & smRs(:,1)<winEndTms(i),2));
                            if winSpd > runThresh
                                actPosn(i) = rad2deg(circ_mean(deg2rad(radPos(smRs(:,1)>=winStartTms(i) & smRs(:,1)<winEndTms(i),2))));
                            end
                        end
                        
                        % Shift everything so laps always show rat going from 0 to 360
                        actPosn = rad2deg(circ_dist(deg2rad(actPosn), deg2rad(rewLoc)));
                        actPosn(actPosn<0) = actPosn(actPosn<0) + 360;
                        
                        
                        % Adjust the number of times of the PPM since it has the potential to be miscalculated by the decoder
                        ppm(:,nWin+1:size(ppm,2)) = [];
                        
                        % Extract the most probable position from PPM
                        decodedPosBins = nan(1,size(ppm,2));
                        decodedPosns = nan(1,size(ppm,2));
                        for i = 1:size(ppm,2)
                            if sum(sum(diff(ppm(:,i)))) ~= 0 %if probability wasn't equal for all bins (no spikes)
                                tmpPpm = ppm(:,i);
                                [~,decodedPosBins(i)] = max(tmpPpm);
                            end
                        end
                        decodedPosns(~isnan(decodedPosBins)) = radBinCtrs(decodedPosBins(~isnan(decodedPosBins)));
                        
                        
                        %Calculate error of decoding for this lap
                        decodingErr = abs(rad2deg(circ_dist(deg2rad(actPosn), deg2rad(decodedPosns))));
                        dayDecodeErr = [dayDecodeErr decodingErr(~isnan(decodingErr))]; %#ok
                        avgErr(lp) = rad2deg(circ_mean(deg2rad(decodingErr(~isnan(decodingErr)))'));
                        semErr(lp) = rad2deg(circ_std(deg2rad(decodingErr(~isnan(decodingErr)))'));
                        
                        if plotActVsDecPosnEachLap == 1
                            % Plot rat's actual vs decoded position for each lap
                            figure(lapFig);
                            if lp < 24
                                subplot(4,6,lp);
                            end
                            tmVctr = winMidTms-winStartTms(1);
                            imagesc(tmVctr, radBinCtrs, ppm);
                            axis xy;
                            colormap jet;
                            hold on;
                            ln = line([tmVctr(1) tmVctr(end)], [180 180]);
                            set(ln, 'Color', [.5 .5 .5], 'LineWidth', 2);
                            plot(tmVctr(~isnan(decodedPosns)), decodedPosns(~isnan(decodedPosns)), '.', 'Color', [1 1 0]);
                            plot(tmVctr, actPosn, 'LineWidth', 2, 'Color', rgb('DarkOrange'));
                            title(['Lap ' num2str(lp)]);
                            xlabel('Time (s)');
                            if mod(lp,6)==1
                                ylabel('Position (degrees)');
                            end
                        end %if plotting
                        
                        
                        %Convert actual positions to the closest index for the next step
                        posInds = match(actPosn,radBinCtrs);
                        posInds(isnan(actPosn)) = NaN;
                        
                        
                        for p = 1:length(radBinCtrs)
                            % Add up the probability the rat is at each location for the time windows in which the rat’s true position is being indexed
                            confMat(:,p) = confMat(:,p) + nansum(ppm(:,posInds==p),2);
                            %The number of time windows for which the rat's true position was equal to the indexed position
                            posCount(p) = posCount(p) + nansum(posInds==p);
                        end
                        %**NOTE: If there isn't a single spike, the location at which the rat spent the most time will have the greatest values in the confusion matrix
                        %          - The decoder will say probability is equal for all locations if no spiking occurred during a time window.
                        %          - If we then add up all those time windows split by where the rat actually was, the spot that he occupied the most
                        %            will have the highest values.
                        
                        
                    end %laps
                    
                    if saveLapPlots == 1  &&   plotActVsDecPosnEachLap == 1
                        figure(lapFig);
                        curDir = pwd;
                        cd(saveDir1);
                        figName = [group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_begin' num2str(b)];
                        print(figName, '-dpng');
                        cd(curDir);
                        close(lapFig);
                    end
                    
                    
                    if plotErrByBegin == 1
                        %Plot error by lap for each 'begin'
                        figure(errFig);
                        subplot(1,4,b);
                        plot(1:length(avgErr), avgErr, 'k');
                        axis square;
                        xlabel('Lap #');
                        fix_font;
                        title(['Begin ' num2str(b)]);
                        set(gca, 'XTick', 1:length(avgErr));
                        xlim([0.5 length(avgErr)+.5]);
                        if b == 1
                            ylabel({'Average Decoding'; 'Error (degrees)'});
                        end
                        yBnds = get(gca, 'YLim');
                        yMaxes(b) = yBnds(2);
                    end
                    
                    
                    
                end %begin
                
                if plotErrByBegin == 1
                    % Make each 'begin' subplot have the same y scale
                    figure(errFig);
                    for b = 1:4
                        subplot(1,4,b);
                        ylim([0 max(yMaxes)]);
                    end
                    
                    if saveBeginErrPlots == 1
                        cd(saveDir2);
                        figName = [group(g).name '_' group(g).rat(r).name '_day' num2str(d)];
                        print(figName, '-dpng');
                        cd(curDir);
                        close(errFig);
                    end
                end
                
                % Plot cumulative proportion of error for each session
                [cumProp, xScale] = calc_cum_prop_ip_range(dayDecodeErr, [0 180]);
                figure(dayDecErrFig);
                plot(xScale, cumProp, 'Color', [.3 .3 .3]);
                
                % Calculate the confusion matrix for this day
                confMat = confMat./repmat(posCount,1, length(radBinCtrs));
                confMat(isinf(confMat)) = 0;
                
                % Concatenate across days
                allConfMats = cat(3, allConfMats, confMat);
                
                
                % Plot and save for individual days
                if plotConfMatsByDay ==1
                    figName = [group(g).name '_' group(g).rat(r).name '_day' num2str(d)];
                    figure('name', figName);
                    colMap = define_cust_color_map('white', 'black', 64);
                    imagesc(radBinCtrs, radBinCtrs, confMat)
                    colormap(colMap);
                    title({['Day ' num2str(d)]; ['#units = ' num2str(size(spkRstr,1))]})
                    axis xy;
                    xlabel('Actual Position (degrees)');
                    ylabel('Decoded Position (degrees)');
                    axis square;
                    colorbar;
                    if saveDayPlots == 1
                        curDir = pwd;
                        cd(saveDir3);
                        print(figName, '-dpng');
                        cd(curDir);
                    end
                end
                
                
            end %if there are more than 2 units
            
        end %day
        
    end %rat
    
    figure(confMatFig);
    subplot(1,2,g);
    colMap = define_cust_color_map('white', 'black', 64);
    imagesc(radBinCtrs, radBinCtrs, nanmean(allConfMats,3))
    colormap(colMap);
    axis xy;
    xlabel('Actual Position (degrees)');
    ylabel('Decoded Position (degrees)');
    title(group(g).name);
    axis square;
    colorbar;
    
end %group


if saveAllSesConfMats == 1
    curDir = pwd;
    cd(saveDir4);
    figName = 'confusionMatrixAcrossSessions_lapMethod';
    print(figName, '-dpng');
    savefig(figName);
    cd(curDir);
end

if saveCumErrPlots == 1
    curDir = pwd;
    cd(saveDir5);
    figure(dayDecErrFig);
    figName = 'cumulativeErrPlot_lapMethod';
    print(figName, '-dpng');
    savefig(figName);
    cd(curDir);
end

end %fnctn






