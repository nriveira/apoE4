function fmr1CircTrack_x_plotBayesianDecoderExample(group)
% function fmr1CircTrack_x_plotBayesianDecoderExample(group)
%
% PURPOSE:
%   Plot examples of Bayesian decoding aross a lap with actual position of
%   the rat.
%
% INPUT:
%   group struct
%
% OUTPUT:
%   Fn-n: Plot showing above for first 8 laps in each begin.
%
% OPTIONS:
%   saveOrNot: save figs (1) or don't (0)
%   plotRewLocs: whether (1) or not (0) to include the reward locations on
%       the plot as dashed lines
%
% MMD
% 7/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 0;
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\DECODING_GENERAL\BayesianDecodingExamples\decodebyLap';
curDir = pwd;

plotRewLocs = 1; %whether or not to plot the reward locations with dashed lines

%% INITIALIZE

rmBinSz = 4; %same as used elsewhere

radBinCtrs = group(2).rat(1).day(1).binCtrs;

newRewLoc = 0; %used to get all plots to line up - lap goes from reward - reward
[~,newRewInd] = min(abs(circ_dist(deg2rad(radBinCtrs), deg2rad(newRewLoc))-0));

sampRate = 20000; %HzsampRate = 20000; %Hz - spike sampling rate
% bayesWin = 40/1000;
% bayesStep = 10/1000;
bayesWin = 250/1000;
bayesStep = 100/1000;

colMap = define_cust_color_map('white', 'black', 30);

cd(saveDir)

%% GET DATA/MAKE FIGS

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            badU = [];
            uIDs = zeros(length(group(g).rat(r).day(d).xAllBeginUnitInfo),2);
            rateMaps = zeros(length(group(g).rat(r).day(d).xAllBeginUnitInfo), 360/rmBinSz);
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                if max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).rateMap)>=1 %unit is bad if max firing rate in bin does not exceed 1
                    rateMaps(u,:) = group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap; %Smoothed ratemap
                    uIDs(u,:) = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                else
                    badU = [badU u]; %#ok
                end
            end
            rateMaps(badU,:) = [];
            rateMaps(rateMaps==0) = 0.0001; %get rid of zeros because our Bayesian decoder can't handle 'em.
            uIDs(badU,:) = [];
            
            for b = 1:4
                figtitle = [group(g).rat(r).name '_D' num2str(d) 'B' num2str(b)];
                figure('Name', figtitle, 'Position', [270 158 1392 799])
                for ll = 1:length(group(g).rat(r).day(d).begin(b).lapTms)
                    if ll <= 8
                        subplot(4,2,ll)
                        lapStart = group(g).rat(r).day(d).begin(b).lapTms(ll,1);
                        lapEnd = group(g).rat(r).day(d).begin(b).lapTms(ll,2);
                        lapDur = lapEnd - lapStart;
                        
                        if group(g).rat(r).day(d).rewLocs(1) ~= 0
                            [~,rewInd] = min(abs(circ_dist(deg2rad(radBinCtrs), deg2rad(group(g).rat(r).day(d).rewLocs(1)))-0));
                            shiftVal = newRewInd - rewInd;
                            
                        end
                        nTimeBins = round(lapDur * sampRate); %number of bins in spike raster
                        spkRstr = zeros(size(uIDs,1), nTimeBins);
                        
                        timeaxis = 0:bayesStep:lapDur;
                        
                        uCntr = 0; %can't just use u since we discard some units due to firing rate (or down sampling)
                        for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                            
                            uID = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                            if ismember(uID, uIDs, 'row') %if the unit wasn't discarded due to low firing ratemap
                                uCntr = uCntr + 1;
                                
                                allSpkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                                lapSpkTms = allSpkTms(allSpkTms>=lapStart & allSpkTms<=lapEnd);
                                
                                % Take spike times and fill in the raster
                                timePassed = lapSpkTms - lapStart;
                                spkInds = round(timePassed * sampRate);
                                spkInds(spkInds==0)=1;
                                
                                spkRstr(uCntr, spkInds) = 1;
                                
                            end %if unit wasn't discarded
                        end %units
       
                        ppm = BayesianDecoder(spkRstr,rateMaps,bayesWin,bayesStep,sampRate); %Ernie's decoder
                        ppm(isnan(ppm)) = 1/size(rateMaps,2); %nan where pxn is chance
                        if group(g).rat(r).day(d).rewLocs(1) ~= 0
                            ppm = circshift(ppm, shiftVal, 1); %shift ppm in the 1st dimension whatever the shift val is
                        end
                        
                        rpStart = find(group(g).rat(r).day(d).begin(b).radPos(:,1) >= lapStart, 1, 'First');
                        rpEnd = find(group(g).rat(r).day(d).begin(b).radPos(:,1) >= lapEnd, 1, 'First');
                        lapRadPos = group(g).rat(r).day(d).begin(b).radPos(rpStart:rpEnd,2);
                        if group(g).rat(r).day(d).rewLocs(1) ~= 0
                            lapRadPos = lapRadPos - (group(g).rat(r).day(d).rewLocs(1) - newRewLoc); %shift this based on reward location too
                            lapRadPos = wrapTo360(lapRadPos); %make those neg values within 0-360 degrees
                        end
                        
                        lapRadPos = lapRadPos / rmBinSz;
                        
                        imagesc(0:lapDur, 0:360/rmBinSz, ppm)
                        axis xy
%                         colormap(parula)
colormap(colMap)
                        caxis([0 0.2])
                        
                        hold on;
                        plot(0:lapDur/length(lapRadPos):lapDur-(lapDur/length(lapRadPos)), lapRadPos, 'Color', rgb('Red'))
                        
                        if plotRewLocs == 1
                            line([0 size(ppm,2)], [size(ppm,1) size(ppm,1)], 'LineStyle', '--', 'Color', 'White') %plot at the "0" reward location
                            
                            if length(group(g).rat(r).day(d).rewLocs) == 2 %some rats only had 1 reward location? Maybe 1 rat for 1 day
                                line([0 size(ppm,2)], [180/rmBinSz 180/rmBinSz], 'LineStyle', '--', 'Color', 'White') %plot second at 180
                            end %if 2 reward locations
                            
                        end %plot reward locations
                        
                        xticks(0:10:lapDur) %labels every 10 seconds
                        xlabel('Time (s)')
                        
                        yticks(0:18:360/rmBinSz)
                        yticklabels({0, 90, 180, 270, 360})
                        ylabel('Position (deg)')
                        
                        ttl = [group(g).rat(r).name ' D' num2str(d) 'B' num2str(b) 'L' num2str(ll)];
                        title(ttl)
                       
                    end %laps 1-10
                end %laps
                
                cbr = colorbar;
                set(cbr, 'Position', [.92 .775 .01 .15])
                ylabel(cbr, 'Probability')
                keyboard
                if saveOrNot == 1
                    saveas(gcf, figtitle, 'epsc');
                    saveas(gcf, figtitle, 'fig');
                    saveas(gcf, figtitle, 'png');
                end
            end %begin
        end %day
        close all
    end %rat
end %group

cd(curDir)

end %function