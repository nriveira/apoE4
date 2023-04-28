function fmr1CircTrack_x_plotNormRateMapsbyRadPos(group)
% function fmr1CircTrack_x_plotNormRateMapsbyRadPos(group)
%
% PURPOSE:
%  Function plots all units' ratemaps by radial position.
%
% INPUT:
%  group = main data structure for the project
%
% OUTPUT:
%   F1-n: Figure for each day with the normalized firing rate of each cell
%       plotted by position. Identified place fields are noted by the colored
%       dots. Reward zone is noted by the grey lines.
%
% JB Trimper - edited MMD
% 12/8/20
% Colgin Lab

% EDIT: Limited number of cells that are plotted in each figure to improve
% legibility. Added columns. MMD 7/2021

%% OPTIONS

saveOrNot = 0; % Set to 1 to save plots in saveDir (1 fig for each group/rat/day)
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\1DRateMapsWithPlaceFieldsMarked';

maxCol = 20; %max number of cells in a column
maxCell = maxCol * 3; %max number of cells in a plot, for legibility

%% INITIALIZE

binSize = 4; %degrees
binCtrs = linspace(binSize/2, 360-binSize/2, 360/binSize);
newRewLoc = 90;
[~,newRewInd] = min(abs(circ_dist(deg2rad(binCtrs), deg2rad(newRewLoc))-0));

curDir = pwd;
cd(saveDir)

%% GET AND PLOT DATA

for g = 1:2
    fprintf('Group %d\n', g);
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            
            % Find how much to shift the ratemaps to put reward locations at 90 & 270
            rewLoc = group(g).rat(r).day(d).rewLocs(1);
            [~,rewInd] = min(abs(circ_dist(deg2rad(binCtrs), deg2rad(rewLoc))));
            shiftVal = newRewInd - rewInd;
            
            numU = length(group(g).rat(r).day(d).xAllBeginUnitInfo);
            numFigs = ceil(numU/maxCell);
            
            uCntr = 0;
            for f = 1:numFigs
                
                if numFigs == 1
                    figtitle = [group(g).name '_' group(g).rat(r).name '_day' num2str(d)];
                else
                    figtitle = [group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_' num2str(f) 'of' (num2str(numFigs))];
                end %how many figs
                
                figU = numU-uCntr; %units for this figure
                if figU > maxCell
                    figU = maxCell - figU;
                end
                
                numCol = ceil(figU/maxCol);
                uPerCol = zeros(1,numCol);
                uPerCol(numCol) = figU - (numCol-1)*maxCol;
                if numCol >1
                    uPerCol(1,1:numCol-1) = maxCol;
                end
                
                figure('Position', [[70 85 574*numCol 901]], 'name', figtitle);
                hold on;
                
                for col = 1:numCol
                    subplot(1,numCol,col)
                    
                    % Mark reward locations with gray lines
                    ln = line([newRewLoc newRewLoc], [.45 numU + 0.45]);
                    set(ln, 'Color', [0.7  0.7  0.7], 'LineWidth', 5);
                    ln = line([newRewLoc+180 newRewLoc+180], [.45 numU + 0.45]);
                    set(ln, 'Color', [0.7  0.7  0.7], 'LineWidth', 5);
                    alpha .3
                    hold on;
                    
                    yVal = 1;
                    for u = 1:uPerCol(col)
                        uCntr = uCntr + 1;
                        % Get the smoothed ratemap, then normalize between 0 & 1
                        origSmRm = group(g).rat(r).day(d).xAllBeginUnitInfo(uCntr).smRateMap;
                        pkFr = max(origSmRm);
                        smRm = origSmRm ./ pkFr;
                        smRm = rescale([0 smRm],yVal-.5, yVal+.4);
                        smRm(1) = [];
                        
                        % Shift the ratemaps so reward locations are always at 90 & 270
                        shiftedSmRm = circshift(smRm, shiftVal);
                        
                        % Plot the rate map (tuning curve)
                        plot(binCtrs, shiftedSmRm);
                        ln = line([binCtrs(1) binCtrs(end)], [yVal-.5 yVal-.5]);
                        set(ln, 'Color', [.4 .4 .4], 'LineStyle', '--');
                        text(binCtrs(end-7), yVal, [num2str(round(pkFr,2)) ' Hz']);
                        
                        % If there are placefields for this rat, mark those on the tuning curve
                        if ~isempty(group(g).rat(r).day(d).xAllBeginUnitInfo(uCntr).pf)
                            for p = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo(uCntr).pf)
                                
                                %Shift the place field location to line up with the shifted rate-map
                                pfInds = group(g).rat(r).day(d).xAllBeginUnitInfo(uCntr).pf(p).inds;
                                tmpVctr = zeros(1,length(smRm));
                                tmpVctr(pfInds) = smRm(pfInds);
                                pfRates = circshift(tmpVctr, shiftVal);
                                
                                plot(binCtrs(pfRates~=0), pfRates(pfRates~=0), '.', 'MarkerSize', 10);
                                
                            end %place fields
                        end %if there are place fields for this unit
                        
                        yVal = yVal + 1;
                    end %unit
                    
                    % Formatting & labelling
                    ylim([.45 maxCol+.45]);
                    xlim([0 360]);
                    yticks(1:uPerCol(col))
                    if uPerCol(col) == maxCol
                        yticklabels(uCntr-19:uCntr)
                    else
                        yticklabels(uCntr-uPerCol(col)+1:uCntr)
                    end %what to label
                    xlabel('Position (degrees)');
                    ylabel('Unit ID');
                    title({'Firing Rate By Position'; 'For Each Unit'});
                    fix_font;
                    
                end %column
                % Save the plots, if desired
                if saveOrNot == 1
                    saveas(gcf, figtitle, 'epsc');
                    saveas(gcf, figtitle, 'fig');
                    saveas(gcf, figtitle, 'png');
                end %save or not
            end %figs
        end %day
    end %rat
end %group

cd(curDir)
end %fnctn