function fmr1CircTrack_x_acrossLap_alignedFiringRate(group)
% function fmr1CircTrack_x_acrossLap_alignedFiringRate(group)
% 
% PURPOSE:
%   Plot aligned firing rate maps for cells across laps for each day. Can
%   be used to visualize place cell stability across the day.
% 
% INPUT:
%   group struct
% 
% OUPUT:
%   F1-n: Figure for each day with the heatmap (white-black) showing
%       normalized firing rate across position (x-axis) by cell (y-axis)
%       for each lap.
% 
% MMD
% 8/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 1;
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\alignedFiringRateMaps_acrossLaps';

%% INITIALIZE

spatBinSz = 4;
velFilt = 1;
durCrit = 1;

cols = {'Blue', 'Red'};

curDir = pwd;
cd(saveDir)

%% GET DATA

cd(saveDir)
for g = 1:2
    cd(group(g).name)
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            figtitle = [group(g).rat(r).name '_' group(g).rat(r).day(d).name];
            figure('Name', figtitle, 'Position', [495 116 854 866])
            
            pkPos = zeros(1,length(group(g).rat(r).day(d).xAllBeginUnitInfo));
            uMaxes = zeros(1,length(group(g).rat(r).day(d).xAllBeginUnitInfo));
            uIDs = zeros(length(group(g).rat(r).day(d).xAllBeginUnitInfo),2);
            badUs = []; %initialize
            
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                [maxVal, maxInd] = max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap);
                pkPos(u) = group(g).rat(r).day(d).binCtrs(maxInd);
                
                uIDs(u,:) = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                if maxVal < 1
                    badUs = [badUs; u];
                end %will remove later
                
                for b = 1:4
                    if max(group(g).rat(r).day(d).begin(b).unit(u).smRateMap) > uMaxes(u)
                        uMaxes(u) = max(group(g).rat(r).day(d).begin(b).unit(u).smRateMap);
                    end
                end %begin
            end %unit
            
            uIDs(badUs,:) = [];
            [~, sortOrd] = sort(pkPos);
            if ~isempty(uIDs)
                for b = 1:4
                    subplot(4,1,b)
                    xlim([0 12])
                    try
                        ylim([0 length(uIDs)])
                    catch
                        keyboard
                    end
                    
                    lapCntr = 0;
                    lapTms = group(g).rat(r).day(d).begin(b).lapTms; %shorten
                    
                    radPos = group(g).rat(r).day(d).begin(b).radPos;
                    coords = group(g).rat(r).day(d).begin(b).coords;
                    
                    for lp = 1:size(lapTms,1)
                        lapCntr = lapCntr + 1;
                        lapVect = zeros(length(pkPos),360/spatBinSz);
                        
                        startPos = find(radPos(:,1) >= lapTms(lp,1), 1, 'First');
                        endPos = find(radPos(:,1) <= lapTms(lp,2), 1, 'Last');
                        
                        cutRadPos = radPos(startPos:endPos,:);
                        curCoords = coords(startPos:endPos,:);
                        
                        for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                            uID = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                            if ismember(uID, uIDs, 'row') %if the unit wasn't discarded due to low firing ratemap
                                spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                                cutSpkTms = spkTms(spkTms >= lapTms(lp,1) & spkTms <= lapTms(lp,2));
                                
                                rateMap = get_ratemap_circtrack(cutSpkTms, curCoords, cutRadPos, spatBinSz, velFilt, durCrit);
                                smRateMap = smooth_circtrack_ratemap(rateMap, spatBinSz);
                                
                                normMap = smRateMap ./ uMaxes(u);
                                lapVect(sortOrd == u,:) = normMap;
                                
                            end %if unit was not discarded
                        end %unit
                        
                        lapVect(badUs,:) = [];
                        xAx = lapCntr - 1:lapCntr;
                        
                        hold on;
                        imagesc(xAx, 0:length(uIDs), lapVect)
                        caxis([0 1])
                        axis xy
                        colormap(flipud(gray))
                    end %lap
                    
                    for lp = 1:size(lapTms,1)
                        line([lp lp], [0 length(uIDs)], 'Color', rgb(cols{g}), 'LineStyle', '--')
                    end %lp
                    
                    xlabel('Lap number')
                    xticks(1:size(lapTms,1))
                    ylabel({['Begin ' num2str(b)], 'Cell number'})
                    
                end %begin
            end %there are any cells
            cbr = colorbar;
            cbr.Position = [ 0.92    0.7761    0.02    0.1489];
            ylabel(cbr, 'Firing rate')
            cbr.TickLabels = {'Min', '', 'Max'};
            
            if saveOrNot == 1
                saveas(gcf, figtitle, 'epsc');
                saveas(gcf, figtitle, 'fig');
                saveas(gcf, figtitle, 'png');
            end %save or not
        end %day
    end %rat
    cd ../
end %group

cd(curDir)
end %function