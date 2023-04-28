function fmr1CircTrack_x_betweenLap_popVectCorr(group)
% function fmr1CircTrack_x_betweenLap_popVectCorr(group)
% 
% PURPOSE:
%   Calulate and plot population vector correlation by lap difference for
%       the WT and KO rats (as in Hwaun & Colgin 2019 - Fig. 2A).
% 
% INPUT:
%   group = data struct.
% 
% OUTPUT:
%   F1: Population vector correlation by lap difference.
% 
% OPTIONS:
%   minCell = minimum number of cells in a day for a day to be included.
%       Note: Fig will save with minCell in title, so multiple can be
%       plotted.
% 
% MMD
% 8/2021
% Colgin Lab

%% OPTIONS

saveOrNot =1;
saveDir = 'E:\Rat381\results\placecell_props';

minCell = 20;

%% INITIALIZE

popVectCorr = cell(1,1); %by group

spatBinSz = 4;
velFilt = 1;
durCrit = 1;

cols = {'Blue'};

curDir = pwd;
cd(saveDir)

%% GET DATA

for g = 1
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            if length(group(g).rat(r).day(d).xAllBeginUnitInfo) > minCell
                acLapMatrix = {}; %initialize
                uMaxes = zeros(1,length(group(g).rat(r).day(d).xAllBeginUnitInfo));
                uIDs = zeros(length(group(g).rat(r).day(d).xAllBeginUnitInfo),2);
                
                for b = 1:4
                    for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                        if max(group(g).rat(r).day(d).begin(b).unit(u).smRateMap) > uMaxes(u)
                            uMaxes(u) = max(group(g).rat(r).day(d).begin(b).unit(u).smRateMap);
                        end
                        if b == 1
                            uIDs(u,:) = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                        end %just to grab unit IDs
                    end %unit
                end %begin
                uIDs(uMaxes < 1,:) = []; %must reach firing rate of at least 1 Hz to be included
                
                lapCntr = 0; %initialize
                for b = 1:4
                    lapTms = group(g).rat(r).day(d).begin(b).lapTms; %shorten
                    
                    radPos = group(g).rat(r).day(d).begin(b).radPos;
                    coords = group(g).rat(r).day(d).begin(b).coords;
                    
                    for lp = 1:size(lapTms,1)
                        
                        lapCntr = lapCntr + 1;
                        lapVect = [];
                        
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
                                lapVect = [lapVect; normMap];
                                
                            end %if unit was not discarded
                        end %unit
                        
                        acLapMatrix{lapCntr} = lapVect;
                        
                        if length(acLapMatrix) > 1
                            for c = 1:length(acLapMatrix) - 1
                                compVect = acLapMatrix{c};
                                
                                scMat = corrcoef(lapVect, compVect);
                                popVectCorr{g} = [popVectCorr{g}; lapCntr-c scMat(2)];
                                
                            end %compare
                        end %mult laps to compare
                    end %laps
                end %begin
            end %at least 20 cells
        end %day
    end %rat
end %group

%% FIGS

figtitle = ['PopulationVectorCorrelation_acrossLaps_minCell' num2str(minCell)];
figure('Name', figtitle, 'Position', [516 544 904 366])

for g = 1
    subplot(1,2,g)
    
    plot(popVectCorr{g}(:,1), popVectCorr{g}(:,2), '.', 'Color', rgb(cols{g}))
    
    maxDiff = max(popVectCorr{g}(:,1));
    
    fitLine = polyfit(popVectCorr{g}(:,1), popVectCorr{g}(:,2), 1);
    xFit = 1:maxDiff;
    yFit = xFit * fitLine(1) + fitLine(2);
    
    hold on;
    plot(xFit, yFit, 'k')
    
    scMat = corrcoef(popVectCorr{g}(:,1), popVectCorr{g}(:,2));
    R2 = scMat(2)^2;
    
    txt = ['slope = ' num2str(round(fitLine(1),2)) ', R^2 = ' num2str(round(R2,2))];
    text(1, 0.1, txt)

    xlabel('Lap difference')
    xlim([0 maxDiff])
    ylabel('Population vector correlation')
    ylim([0 1])
    
    title(group(g).name)
end %group

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

cd(curDir)
end %function