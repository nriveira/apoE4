function adCircTrack_x_plotTheta(group)
% function fmr1CircTrack_x_plotTheta(group)
%
% PURPOSE:
%  Function to plot LFPs from each tetrode for a single lap for the sake of comparing across
%  to identify which tetrode will serve as the phase reference for future analyses.
%
% INPUT:
%   group = project main data struct
%
% OUTPUT:
%   F1-N: plots, one for each tetrode showing raw LFP, narrow theta 
%       (6-10 Hz) filtered, and broad (3-20 Hz) filtered
% 
% NOTE MMD: This project has changed quite a bit since John made did the
%   preliminary analyses, and the theta tet isn't as wildly used in its
%   new iteration. 
%
% JB Trimper
% 2/23/21
% Colgin Lab

%% INITIALIZE

curDir = pwd; %Note where we are so we can return to it later
b = 1; %Begin 1 only
lpNum = 2; %Lap to plot the data from (arbitrary; just an option to look at data from dif laps)
indsToPlot = [1 6000]; %which data from the lap, cuz looking at a whole lap's data is too much for the small axis

saveOrNot = 0; %just in case I need to show Laura or something

saveDir = 'E:\Rat381\results\plotTheta';

%% GET DATA/MAKE PLOTS

for g = 1
    fprintf('Group %d\n', g);
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            
            uIDs = [];
            for u = 1:length(group(g).rat(r).day(d).begin(b).unit)
                uIDs(u) = group(g).rat(r).day(d).begin(b).unit(u).ID(1); %#ok
            end
            
            tetNums = unique(uIDs, 'stable');
            
            uCounts = histcounts(uIDs);
            uCounts(uCounts==0)=[];
            
            %Start and end time for this 'begin'
            startLapTm = group(g).rat(r).day(d).begin(b).lapTms(lpNum,1);
            endLapTm =  group(g).rat(r).day(d).begin(b).lapTms(lpNum,2);
            
            cd(group(g).rat(r).day(d).begin(b).dir)
            figName = [group(g).name '_' group(g).rat(r).name '_Day' num2str(d)];
            figure('name', figName, 'Position', [1 31 1920 973]);
            
            yBnds = [];
            for tt = 1:length(tetNums)
                
                lfpRoot = ['CSC' num2str(tetNums(tt))];
                lfpStruct = read_in_lfp(['CSC' num2str(tetNums(tt)) '.ncs']);
                
                load([lfpRoot '_broadThetaLfp.mat']); %#ok
                lfpStruct.broadThetaLfp = filtLfp;
                
                load([lfpRoot '_narrowThetaLfp.mat']); %#ok
                lfpStruct.narrowThetaLfp = filtLfp;
                
                lapLfpInds = find(lfpStruct.ts>=startLapTm & lfpStruct.ts<=endLapTm);
                
                subplot(4,3,tt);
                plot(lfpStruct.data(lapLfpInds));
                hold on;
                plot(lfpStruct.broadThetaLfp(lapLfpInds));
                plot(lfpStruct.narrowThetaLfp(lapLfpInds), 'LineWidth', 1.5);
                title({['Tet #' num2str(tetNums(tt))]; [num2str(uCounts(tt)) ' Units']});
                zero_line;
                
                yBnds = [yBnds; get(gca, 'YLim')]; %#ok
                
                xlim(indsToPlot);
                
            end %tetrode
            
            if saveOrNot == 1
                cd(saveDir)
                saveas(gcf, figName, 'epsc')
                saveas(gcf, figName, 'png')
                saveas(gcf, figName, 'fig')
            end %save option
            
        end %day
        
    end %rat
    
end %group

cd(curDir);

end %fnctn