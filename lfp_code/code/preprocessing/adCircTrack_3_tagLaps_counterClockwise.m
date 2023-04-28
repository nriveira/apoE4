function group = adCircTrack_3_tagLaps_counterClockwise(group)
% function group = fmr1CircTrack_3_tagLaps(group)
%
% PURPOSE
%  Functions add indices to main data structure that point to times which the rat was moving from
%  reward location 1 to location 2 (and vice versa) as well as times which the rat was completing
%  a full lap from/to location 1
%
% INPUT:
%  group = main data structure for the project
%
% OUTPUT:
%  A subfield within each 'begin' sub-structure that contains the indices for the start and stop of each
%  outbound trajectory, inbound trajectory, and full laps
%
% JB Trimper
% 12/15/20
% Colgin Lab

plotCheck = 0; %Way to doublecheck that paths are covering the expected portions of the track
%              Makes 2 figures for each begin window (for each day/rat/group.
savePlots = 0;
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\lapSortingPlots';

rewLocInBnd = 15; % +/- rewLocInBnd degrees -- around reward location --> This creates space around the reward location for you to say a lap has started or ended
rewLocOutBnd = 30; % +/- rewLocOutBnd degrees --  around reward location --> This helps us ensure the start/stop were still close to the reward loc (between rewLocInBnd & rewLocOutBnd)

for g = 1:2
    fprintf('Group %d\n', g);
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            
            rewLocs = group(g).rat(r).day(d).rewLocs;
            
            if plotCheck == 1
                if length(rewLocs)>1
                    figName = [group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_HalfLaps'];
                    halfLapFig = figure('name', figName, 'Position', [256 254 1052 738]); %Laptop
                    %             halfLapFig = figure('name', figName, 'Position', [2    32   958   964]); %Desktop
                end
                
                figName = [group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_FullLaps'];
                fullLapFig = figure('name', figName, 'Position', [ 258         366        1051         459]); %Laptop
                %             fullLapFig = figure('name', figName, 'Position', [ 1250 32 670 964]); %Desktop
            end %plot check make plots
            for b = 1:4
                fprintf('\t\t\tBegin %d\n', b);
                
                radPos = group(g).rat(r).day(d).begin(b).radPos;
                coords = group(g).rat(r).day(d).begin(b).coords;
                
                
                %% GET TRAVERSALS BETWEEN REWARD LOCATIONS  ~~~~~> ALWAYS COUNTERCLOCKWISE!
                if length(rewLocs)>1
                    radPosBnry = zeros(1,size(radPos,1));
                    
                    %   Outbound = from location 1 to location 2
                    %     Pt 1: Define trajectory bounds
                    obTrajBnds = [rewLocs(1)+rewLocInBnd rewLocs(2)-rewLocInBnd]; %ob trajectories will go from the + side of rew loc 1 to - side of rew loc 2
                    %     Pt 2: Label times when rat was within bounds with boolean
                    radPosBnry(radPos(:,2)>=obTrajBnds(1) & radPos(:,2)<=obTrajBnds(2)) = 1;
                    
                    %   Inbound = from location 2 to location 1
                    %     Pt 1a: Define trajectory bounds
                    ibTrajBnds = [rewLocs(2)+rewLocInBnd rewLocs(1)-rewLocInBnd]; %ib trajectories will go from the - side of rew loc 2 to + side of rew loc 1
                    %     Pt 1b: Ib trajectories can cross 0/360, so the bounds calculated above need to be adjusted to line up with our data (0-360)
                    if ibTrajBnds(1)>360 && ibTrajBnds(2)>0
                        ibTrajBnds(1) = ibTrajBnds(1) - 360;
                    elseif ibTrajBnds(1)>0 && ibTrajBnds(2)<0
                        ibTrajBnds(2) = ibTrajBnds(2) + 360;
                    end
                    
                    %     Pt 2: Label times within bounds as with bool vals, but account for fluctuation around 0/360
                    if ibTrajBnds(1) < ibTrajBnds(2)
                        radPosBnry(radPos(:,2)>=ibTrajBnds(1) & radPos(:,2)<=ibTrajBnds(2)) = 1;
                    elseif ibTrajBnds(1) > ibTrajBnds(2)
                        radPosBnry(radPos(:,2)>=ibTrajBnds(1) | radPos(:,2)<=ibTrajBnds(2)) = 1;
                    end
                    
                    
                    if plotCheck == 1
                        figure(halfLapFig);
                        spTtls = {'Outbound', 'Inbound', 'Bad Paths'};
                        for i = 1:3
                            subplot(3,4,(i-1)*4+b);
                            hold on;
                            axis([-54 54 -54 54]);
                            axis square;
                            if i == 1
                                title(['Begin ' num2str(b)]);
                            elseif i == 3
                                xlabel('X Position (cm)');
                            end
                            if b == 1
                                ylabel({['Begin ' num2str(b)], 'Y Position (cm)'});
                                ylabel({spTtls{i}, 'Y Position (cm)'});
                            end
                            plot(coords(:,2), coords(:,3), 'Color', [.4 .4 .4]);
                        end
                    end
                    
                    
                    % Sort the boolean vector into actual ob vs ib (vs artefact) trajectories
                    posChunks = bwconncomp(radPosBnry,4); % find consecutive chunks of 1s
                    outCntr = 0; %for counting ib/ob laps
                    inCntr = 0;
                    for c = 1:length(posChunks.PixelIdxList)
                        chunkInds = posChunks.PixelIdxList{c};
                        if length(posChunks.PixelIdxList{c}) > 10
                            startChunkInd = chunkInds(1);
                            endChunkInd = chunkInds(end);
                            
                            % Outbound Trajectories (from rew loc 1 to rew loc 2)
                            if rad2deg(circ_dist(deg2rad(radPos(startChunkInd,2)), deg2rad(rewLocs(1)+rewLocOutBnd))) < 0  && ... %If start position was between rewLocInBnd & rewLocOutBnd degrees above rew loc 1
                                    rad2deg(circ_dist(deg2rad(radPos(endChunkInd,2)), deg2rad(rewLocs(2)-rewLocOutBnd))) > 0 %... and end position was between rewLocInBnd & rewLocOutBnd degrees below rew loc 2
                                
                                outCntr = outCntr + 1; %add the ob trajectory's indices and radial coordinates to the structure
                                group(g).rat(r).day(d).begin(b).outPathInds(outCntr,:) = [startChunkInd endChunkInd];
                                group(g).rat(r).day(d).begin(b).outPathTms(outCntr,:) = [radPos(startChunkInd,1) radPos(endChunkInd,1)];
                                
                                if plotCheck == 1
                                    subplot(3,4,b);
                                    plot(coords(startChunkInd:endChunkInd,2), coords(startChunkInd:endChunkInd,3));
                                end
                                
                                % Inbound Trajectories (from rew loc 2 to rew loc 1)
                            elseif rad2deg(circ_dist(deg2rad(radPos(startChunkInd,2)), deg2rad(rewLocs(2)+rewLocOutBnd))) < 0  && ... %If start position was between rewLocInBnd & rewLocOutBnd degrees above rew loc 2
                                    rad2deg(circ_dist(deg2rad(radPos(endChunkInd,2)), deg2rad(rewLocs(1)-rewLocOutBnd))) > 0 %... and end position was between rewLocInBnd & rewLocOutBnd degrees below rew loc 1
                                
                                inCntr = inCntr + 1; %add the ib trajectory's indices and radial coordinates to the structure
                                group(g).rat(r).day(d).begin(b).inPathInds(inCntr,:) = [startChunkInd endChunkInd];
                                group(g).rat(r).day(d).begin(b).inPathTms(inCntr,:) = [radPos(startChunkInd,1) radPos(endChunkInd,1)];
                                
                                if plotCheck == 1
                                    subplot(3,4,4+b);
                                    plot(coords(startChunkInd:endChunkInd,2), coords(startChunkInd:endChunkInd,3));
                                end
                                
                            else
                                
                                %If the tagged coordinate chunk didn't start near one reward location then end near the other, that chunk isn't a full IB or OB trajectory
                                %  (maybe the rat backed up a few steps or leaned over the bounds then swiveled back to reward)
                                %    We're not saving these 'bad' trajectories, just plotting them if desired
                                if plotCheck == 1
                                    subplot(3,4,8+b);
                                    plot(coords(startChunkInd:endChunkInd,2), coords(startChunkInd:endChunkInd,3));
                                end
                                
                            end
                        end
                    end
                end %if two reward locs
                
                
                
                %% GET FULL LAPS
                %   From reward location 1 all the way around back to reward location 1
                
                innerBnds = [rewLocs(1)+rewLocInBnd  rewLocs(1)-rewLocInBnd]; % [ r1+  (sb1)     r1-  (eb1) ]  <-- these codes are to line up with my notebook
                outerBnds = [rewLocs(1)+rewLocOutBnd  rewLocs(1)-rewLocOutBnd]; % [ r1++ (sb2)     r1-- (eb2) ]
                
                % Based on where the reward location was, the boundaries around them could be <0 or >360. We have to correct for that.
                innerBnds(innerBnds>360) = innerBnds(innerBnds>360)-360;
                innerBnds(innerBnds<0) = innerBnds(innerBnds<0)+360;
                outerBnds(outerBnds>360) = outerBnds(innerBnds>360)-360;
                outerBnds(outerBnds<0) = outerBnds(outerBnds<0)+360;
                
                %Label times where the rat was not at reward 1
                radPosBnry = zeros(1,size(radPos,1));
                if innerBnds(1)>innerBnds(2) %Labeling procedure will vary based on note above the previous chunk of code(Based on where...)
                    radPosBnry(radPos(:,2)>innerBnds(1) | radPos(:,2)<innerBnds(2)) = 1;
                else
                    radPosBnry(radPos(:,2)>innerBnds(1) & radPos(:,2)<innerBnds(2)) = 1;
                end
                
                % From the binary, identify chunks of radial coordinates where the rat was not at reward locations
                lapChnks = bwconncomp(radPosBnry, 4);
                
                if plotCheck == 1
                    %Plot all coordinates on which to superimpose the laps later, for visual inspection
                    figure(fullLapFig);
                    spTtls = {'Full Laps', 'Bad'};
                    for i = 1:2
                        subplot(2,4,(i-1)*4+b);
                        hold on;
                        axis([-54 54 -54 54]);
                        axis square;
                        plot(coords(:,2), coords(:,3), 'Color', [.4 .4 .4]);
                        if i == 1
                            title(['Begin ' num2str(b)]);
                        else
                            xlabel('X Position (cm)');
                        end
                        if b == 1
                            ylabel({['Begin ' num2str(b)], 'Y Position (cm)'});
                            ylabel({spTtls{i}, 'Y Position (cm)'});
                        end
                    end
                end
                
                
                % Get lap info and save it
                lapCntr = 0;
                for c = 1:length(lapChnks.PixelIdxList)
                    lapInds = lapChnks.PixelIdxList{c};
                    
                    % Don't bother asking if it's a lap if there are less than 20 indices (this would mean a full lap lasted <0.66 sec)
                    if length(lapInds) > 20
                        
                        
                        startChunkInd = lapInds(1);
                        endChunkInd = lapInds(end);
                        
                        startLoc = radPos(startChunkInd,2); %start position for this chunk of radial coordinates
                        endLoc = radPos(endChunkInd,2); %end position
                        
                        %If the start location was between the inner & outer bounds...
                        % and the end location was between the inner and outer bounds on the other side of reward 1
                        if rad2deg(circ_dist(deg2rad(outerBnds(1)),deg2rad(startLoc))) <= rewLocInBnd && ...
                                rad2deg(circ_dist(deg2rad(outerBnds(1)),deg2rad(startLoc))) > 0 && ...
                                rad2deg(circ_dist(deg2rad(endLoc), deg2rad(outerBnds(2)))) <= rewLocInBnd && ...
                                rad2deg(circ_dist(deg2rad(endLoc), deg2rad(outerBnds(2)))) > 0
                            
                            
                            lapCntr = lapCntr + 1;
                            group(g).rat(r).day(d).begin(b).lapInds(lapCntr,:) = [startChunkInd endChunkInd];
                            group(g).rat(r).day(d).begin(b).lapTms(lapCntr,:) = [radPos(startChunkInd,1) radPos(endChunkInd,1)];
                            
                            if plotCheck == 1
                                subplot(2,4,b);
                                plot(coords(startChunkInd:endChunkInd,2), coords(startChunkInd:endChunkInd,3));
                            end
                            
                        else
                            if plotCheck == 1
                                subplot(2,4,4+b);
                                plot(coords(startChunkInd:endChunkInd,2), coords(startChunkInd:endChunkInd,3));
                            end
                        end
                        
                    end
                end
                
                
                
            end %begin
            
            if plotCheck == 1 && savePlots == 1
                curDir = pwd;
                cd(saveDir);
                
                if length(rewLocs)>1
                    figure(halfLapFig)
                    print([group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_HalfLaps'], '-dpng');
                    savefig([group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_HalfLaps']);
                end
                
                figure(fullLapFig)
                print([group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_FullLaps'], '-dpng');
                savefig([group(g).name '_' group(g).rat(r).name '_day' num2str(d) '_FullLaps'])
                
                cd(curDir);
            end
            
        end %day
    end %rat
end %group


end %fnctn




