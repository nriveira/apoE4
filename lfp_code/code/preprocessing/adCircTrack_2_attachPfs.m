function group = adCircTrack_2_attachPfs(group)
% function group = fmr1CircTrack_2_attachPfs(group)
%
% PURPOSE:
%  To attach place-field info for each unit to main data structure.
%
% NOTE:
%  group = output of fmr1CircTrack_1_buildDataStruct
%
% OUTPUT:
%  same as input with pf field within each unit's sub-structure
%
% JBT
% 12/3/2020
% Colgin Lab

spatBinSz = 4; %spatial bin size in degrees
velFilt = 1; %velocity filter the spikes while looking for PFs
durCrit = 1; %enforce duration criteria while looking for PFs
minPkFr = 1; %Hz

gWinStd = 8; %degrees, as in Zheng et al. 2021
gWinSz = gWinStd * 2;

for g = 1:2
%for g = 1 %for only apoE4 rat
    fprintf('Group %d\n', g);
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            
            xBeginSpkCnts = zeros(360/spatBinSz, 4, length(group(g).rat(r).day(d).begin(1).unit));
            xBeginTpb = zeros(360/spatBinSz, 4);
            
            xBeginUIDs = [];
            
            for b = 1:4
                fprintf('\t\t\tBegin %d\n', b);
                
                radPos = group(g).rat(r).day(d).begin(b).radPos;
                coords = group(g).rat(r).day(d).begin(b).coords;
                
                for u = 1:length(group(g).rat(r).day(d).begin(b).unit)
                    fprintf('\t\t\t\tUnit %d/%d\n', u, length(group(g).rat(r).day(d).begin(b).unit));
                    spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                    [rateMap, binCtrs, tpb, spkCnts] = get_ratemap_circtrack(spkTms, coords, radPos, spatBinSz, velFilt, durCrit);
                    
                    smRateMap = smooth_circtrack_ratemap(rateMap, spatBinSz, gWinSz, gWinStd);
                    
                    group(g).rat(r).day(d).begin(b).unit(u).rateMap = rateMap;
                    group(g).rat(r).day(d).begin(b).unit(u).smRateMap = smRateMap;
                    group(g).rat(r).day(d).begin(b).unit(u).pf = get_circtrack_pfs(smRateMap, spatBinSz); 
                    group(g).rat(r).day(d).begin(b).unit(u).pf_backwardExp = get_circtrack_pfs_backwardExpansion(smRateMap); 
                    group(g).rat(r).day(d).begin(b).unit(u).spkCnts = spkCnts;
                    
                    xBeginSpkCnts(:,b,u) = spkCnts;
                    xBeginUIDs = [xBeginUIDs; group(g).rat(r).day(d).begin(b).unit(u).ID]; %#ok

                    
                end %unit
                
                group(g).rat(r).day(d).begin(b).tpb = tpb;
                xBeginTpb(:,b) = tpb;
                
            end %begin
            
            allUSpkCntsXBegins = squeeze(sum(xBeginSpkCnts,2));
            allUSpkCntsX3Begins = squeeze(sum(xBeginSpkCnts(:,1:3,:),2));
            xAllBeginTpb = squeeze(sum(xBeginTpb,2));
            x3BeginTpb = squeeze(sum(xBeginTpb(:,1:3),2));
            xAllBeginUIDs = unique(xBeginUIDs, 'row', 'stable');
            group(g).rat(r).day(d).xAllBeginTpb = xAllBeginTpb';
            group(g).rat(r).day(d).x3BeginTpb = x3BeginTpb';
            group(g).rat(r).day(d).binCtrs = binCtrs;
            
            for u = 1:size(allUSpkCntsXBegins,2)
                rateMap = (allUSpkCntsXBegins(:,u) ./ xAllBeginTpb)';
                smRateMap = smooth_circtrack_ratemap(rateMap, spatBinSz, gWinSz, gWinStd);

                % pf = get_circtrack_pfs_m2(rateMap);
                
                pf = get_circtrack_pfs(smRateMap, spatBinSz);
                pf_b = get_circtrack_pfs_backwardExpansion(smRateMap);
                
                badFs = [];
                for p = 1:length(pf_b)
                    if pf_b(p).pkFr < minPkFr
                        badFs = [badFs; p];
                    end %check
                end %place field
                pf_b(badFs) = [];

                badFs = [];
                for p = 1:length(pf)
                    if pf(p).pkFr < minPkFr
                        badFs = [badFs; p];
                    end %check
                end %place field
                pf(badFs) = [];
          
                group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID = xBeginUIDs(u,:); 
                group(g).rat(r).day(d).xAllBeginUnitInfo(u).spkCnts = allUSpkCntsXBegins(:,u)';
                group(g).rat(r).day(d).xAllBeginUnitInfo(u).rateMap = rateMap;
                group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap = smRateMap;
                group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf = pf;
                group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf_backwardExp = pf_b;
                
                rateMap = (allUSpkCntsX3Begins(:,u) ./ x3BeginTpb)';
                smRateMap = smooth_circtrack_ratemap(rateMap, spatBinSz, gWinSz, gWinStd);
                pf = get_circtrack_pfs(smRateMap, spatBinSz);
                pf_b = get_circtrack_pfs_backwardExpansion(smRateMap, spatBinSz);
                
                group(g).rat(r).day(d).x3BeginUnitInfo(u).ID = xBeginUIDs(u,:); 
                group(g).rat(r).day(d).x3BeginUnitInfo(u).spkCnts = allUSpkCntsX3Begins(:,u)';
                group(g).rat(r).day(d).x3BeginUnitInfo(u).rateMap = rateMap;
                group(g).rat(r).day(d).x3BeginUnitInfo(u).smRateMap = smRateMap;
                group(g).rat(r).day(d).x3BeginUnitInfo(u).pf = pf;
                group(g).rat(r).day(d).x3BeginUnitInfo(u).pf_backwardExp = pf_b;

            end %unit
            
        end %day
    end %rat
end %group



end %fnctn