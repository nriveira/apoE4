function group = adCircTrack_5_addSleepInfoToStruct(group)
% function group = fmr1CircTrack_5_addSleepInfoToStruct(group)
%
% PURPOSE:
%   To add the sleep information to the struct.
%
% INPUT:
%   group = data struct, through function fmr1CircTrack_4_...
%
% OUTPUT:
%   group = data struct, with sleep info
%
% MMD
% Colgin Lab

%% INITIALIZE

restPotDim = [12 12]; %10 cm diameter, with room for overhanging rat

dataDir = 'E:\Rat381';

spdCheck = 1;

stdCut = 5; %standard deviations above the mean for rip power cut
durThresh = 0.5; %minimum dur for ripples

%% GET INFO

for g = 1
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            cd([dataDir '\' group(g).rat(r).day(d).name]);
            
            fid = fopen('TTList.txt');
            tmp = textscan(fid, '%s', 'delimiter', '\n');
            unitList = tmp{1};
            fclose(fid);
            
            slpDirInfo = dir('sleep*');
            
            if isempty(slpDirInfo)
                for s = 1:5
                    group(g).rat(r).day(d).sleep(s).coords = [];
                    group(g).rat(r).day(d).sleep(s).unit = [];
                    group(g).rat(r).day(d).sleep(s).ripInds = [];
                    group(g).rat(r).day(d).sleep(s).ripTms = [];
                    
                end %sleeps
            else
                slpFoldNums = zeros(1,length(slpDirInfo));
                for sInd = 1:length(slpDirInfo)
                    foldName = slpDirInfo(sInd).name;
                    s = str2num(foldName(6));
                    slpFoldNums(sInd) = s;
                end
                
                for s = 1:5
                    if find(slpFoldNums == s)
                        
                        fprintf('\t\t\tSleep %d\n', s);
                        group(g).rat(r).day(d).sleep(s).dir = [dataDir '\' group(g).rat(r).day(d).name '\sleep' num2str(s)];
                        cd(group(g).rat(r).day(d).sleep(s).dir)
                        
                        coords = read_in_coords('VT1.nvt', restPotDim(1), restPotDim(2));
                        
                        group(g).rat(r).day(d).sleep(s).coords = coords;
                        
                        tetNums = [];
                        %                         uCntr = 0;
                        for u = 1:length(unitList)
                            try
                                spkTms = Readtfile(unitList{u});
                                %                         goodRead = 1;
                                %                         uCntr = uCntr + 1;
                            catch
                                %                         goodRead = 0;
                                if ~isfile(unitList{u})
                                    spkTms = [];
                                end
                            end
                            
%                             if goodRead == 1
                                spkTms = spkTms ./ 10^4;
                                
                                hyphInd = strfind(unitList{u}, '_');
                                tetNum = str2double(unitList{u}(3:hyphInd-1));
                                extInd = strfind(unitList{u}, '.t');
                                clustNum = str2double(unitList{u}(hyphInd+1:extInd-1));
                                
%                                 group(g).rat(r).day(d).sleep(s).unit(uCntr).ID = [tetNum clustNum];
%                                 group(g).rat(r).day(d).sleep(s).unit(uCntr).spkTms = spkTms;
                                group(g).rat(r).day(d).sleep(s).unit(u).ID = [tetNum clustNum];
                                group(g).rat(r).day(d).sleep(s).unit(u).spkTms = spkTms;
                                
                                tetNums = [tetNums tetNum]; %#ok
                                
%                             end %try to read in the units
                        end %units
                        fprintf('\t\t\t\t%d Units Attached\n', length(unitList));
                    else
                        group(g).rat(r).day(d).sleep(s).coords = [];
                        group(g).rat(r).day(d).sleep(s).unit = [];
                    end %this s has a folder
                end %sleeps
            end %there are sleep folders
        end %day
    end %rat
    
end %group


end %function