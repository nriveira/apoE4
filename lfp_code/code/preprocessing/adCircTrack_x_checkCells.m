function adCircTrack_x_checkCells(group)
% function fmr1CircTrack_x_checkCells(group)
%
% PURPOSE:
%   Checks that no cells have been included that meet the definition of an
%   interneuron (5 Hz) and checks that all cells have a refractory period.
%
% INPUT:
%   group = data struct
% 
% OUPUT:
%   If applicable, code will specify cells that fail either the interneuron
%   or refractory period check.
% 
% MMD
% 11/2021
% Colgin Lab

%% INITIALIZE

intFR = 5; %Hz

%% CHECK DATA

for g = 1
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                checkInt = 0; %initialize as pass
                checkRP = 0; %initialize as pass
                
                exInfo = {};
                for b = 1:4
                    spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                    diffTms = diff(spkTms);
                    if min(diffTms) < 0.0005 %within a ms
                        checkRP = 1;
                        exInfo{length(exInfo)+1} = ['B' num2str(b)];
                    end %RP
                    
                    totTime =  group(g).rat(r).day(d).begin(b).coords(end,1) - group(g).rat(r).day(d).begin(b).coords(1,1);
                    tmpFR = length(spkTms) / totTime;
                    if tmpFR >= intFR
                        checkInt = 1;
                        exInfo{length(exInfo)+1} = ['B' num2str(b)];
                    end %int
                end %begin
                
                for s = 1:5
                    if ~isempty(group(g).rat(r).day(d).sleep(s).unit)
                        spkTms = group(g).rat(r).day(d).sleep(s).unit(u).spkTms;
                        diffTms = diff(spkTms);
                        if min(diffTms) < 0.0005 %within a ms
                            checkRP = 1;
                            exInfo{length(exInfo)+1} = ['S' num2str(s)];
                        end %RP
                        
                        totTime =  group(g).rat(r).day(d).sleep(s).coords(end,1) - group(g).rat(r).day(d).sleep(s).coords(1,1);
                        tmpFR = length(spkTms) / totTime;
                        if tmpFR >= intFR
                            checkInt = 1;
                            exInfo{length(exInfo)+1} = ['S' num2str(s)];
                        end %int
                    end %if there are spikes
                end %sleep
                
                if checkInt == 1
                    fprintf('%s Rat %d Day %d: TT%d_%d failed interneuron check\n', group(g).name, r, d, group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID(1), group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID(2))
%                     for i = 1:length(exInfo)
%                         fprintf('\t%s\n', exInfo{i})
%                     end %extra info
                end %fprintf int
                if checkRP == 1
                    fprintf('%s Rat %d Day %d: TT%d_%d failed RP check\n', group(g).name, r, d, group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID(1), group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID(2))
%                     for i = 1:length(exInfo)
%                         fprintf('\t%s\n', exInfo{i})
%                     end %extra info
                end %fprintf int
            end %unit
        end %day
    end %rat
end %group

end %function