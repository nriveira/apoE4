function fmr1CircTrack_x_lapRasterPlots(group)
% function fmr1CircTrack_x_lapRasterPlots(group)
%
% PURPOSE:
%   Plot example spike raster plots for successive laps across all four
%   begins. Makes plots for all units. Saves in folders for each day for
%   each rat. Only includes spike times while the rat is running > 5 cm/s.
%
% INPUT:
%   group struct
%
% OUTPUT:
%   Fn-n: Figure for each unit with tick marks indicated the rat's position
%       when the cell fired. Each row represents a different lap.
%
% NOTE:
%   Function closes figures internally, so if you want to see them comment
%       out line 118 or add a keyboard. Makes an individual plot for each
%       unit, so be aware it could be a lot of figures.
%
% MMD
% 7/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 1;
% saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\lapRasterPlots';
saveDir = 'E:\Rat381\results\lapRasterPlots';

cols = {'Blue'};

%% INITIALIZE

binSize = 4; %degrees
binCtrs = linspace(binSize/2, 360-binSize/2, 360/binSize);

runThresh = 5; %cm/s
greycol = rgb('LightGray');

cd(saveDir)

%% GET DATA/MAKE PLOTS

for g = 1
    try
        cd(group(g).name)
    catch
        mkdir(group(g).name)
        cd(group(g).name)
    end
    for r = 1:length(group(g).rat)
        if ~isfolder(group(g).rat(r).name)
            mkdir(group(g).rat(r).name)
        end %if is folder
        cd(group(g).rat(r).name)
        for d = 1:length(group(g).rat(r).day)
            
            if ~isfolder(group(g).rat(r).day(d).name)
                mkdir(group(g).rat(r).day(d).name)
            end %if is folder
            
            cd(group(g).rat(r).day(d).name)
            
            rewLoc = group(g).rat(r).day(d).rewLocs(1);
            shiftDeg = 0 - rewLoc;
            shiftVal = round(shiftDeg/binSize); % the divisor is the radial bin size in degrees

            
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                tmpID = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                figtitle = ['TT' num2str(tmpID(1)) '_' num2str(tmpID(2))];
                
                figure('Name', figtitle)
                
                yVal = 1;
                lpCntr = 0;
                
                for b = 1:4
                    
                    spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                    radPos = group(g).rat(r).day(d).begin(b).radPos;
                    coords = group(g).rat(r).day(d).begin(b).coords;
                    
                    instRs = get_runspeed(coords);
                    smRs = smooth_runspeed(instRs);
                    
                    for lp = 1:size(group(g).rat(r).day(d).begin(b).lapTms,1)
                        lapStart = group(g).rat(r).day(d).begin(b).lapTms(lp,1);
                        lapEnd = group(g).rat(r).day(d).begin(b).lapTms(lp,2);
                        lapSpks = spkTms(spkTms >= lapStart & spkTms <= lapEnd);
                        
                        %                         lapSpkPos = [];
                        
                        for st = 1:length(lapSpks)
                            spkTm = lapSpks(st);
                            posInd = match(spkTm, smRs(:,1));
                            
                            if smRs(posInd,2) < runThresh
                                continue %to next spike
                            end %doesn't reach run thresh
                            
                            %                             lapSpkPos = [lapSpkPos radPos(posInd,2)];
                            shiftedPos = wrapTo360(radPos(posInd,2)+shiftDeg);
                            line([shiftedPos shiftedPos], [yVal-.3 yVal+.3], 'Color', rgb(cols{g}))
                        end %spktms
                        
                        yVal = yVal + 1;
                    end %lap
                    
                    ln = line([binCtrs(1) binCtrs(end)], [yVal-.5 yVal-.5]);
                    set(ln, 'Color', [.4 .4 .4], 'LineStyle', '--');
                    
                    text(binCtrs(end)+5, (size(group(g).rat(r).day(d).begin(b).lapTms,1)/2)+lpCntr, ['Begin ' num2str(b)]);
                    lpCntr = lpCntr + size(group(g).rat(r).day(d).begin(b).lapTms,1);
                end %begin
                
                xlim([0 360])
                xlabel('Position relative to reward 1 (deg)')
                
                ylim([0.45 yVal+0.45])
                ylabel('Lap Number')
                
                title(['Unit: ' 'TT' num2str(tmpID(1)) '\_' num2str(tmpID(2))])
                %                 pf = group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf;
                
                shiftRm = circshift(group(g).rat(r).day(d).xAllBeginUnitInfo(u).rateMap, shiftVal);
                shiftSmRm = circshift(group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap, shiftVal);
                pf = get_circtrack_pfs_m2(shiftRm, shiftSmRm);
                
                for p = 1:length(pf)
                    if pf(p).radPos(1) < pf(p).radPos(end)
                        x = [pf(p).radPos(1) pf(p).radPos(1) pf(p).radPos(end) pf(p).radPos(end)];
                        y = [0 yVal+0.45 yVal+0.45 0];
                        createpatch = patch(x, y, [greycol]);
                        alpha(0.5);
                    else %need to account for place fields that wrap around the 0-360 degree mark
                        x = [0 0 pf(p).radPos(end) pf(p).radPos(end)];
                        y = [0 yVal+0.45 yVal+0.45 0];
                        createpatch = patch(x, y, [greycol]);
                        alpha(0.1);
                        hold on;
                        x = [pf(p).radPos(1) pf(p).radPos(1) 360 360];
                        y = [0 yVal+0.45 yVal+0.45 0];
                        createpatch = patch(x, y, [greycol]);
                        alpha(0.1);
                    end %if
                end %place field
                
                if saveOrNot == 1
                    saveas(gcf, figtitle, 'epsc')
                    saveas(gcf, figtitle, 'png')
                    saveas(gcf, figtitle, 'fig')
                end %save option
                
            end %unit
            close all
            cd ../
            
        end %day
        cd ../
    end %rat
    cd ../
end %group

cd(curDir)
end %function