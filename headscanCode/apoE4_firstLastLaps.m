function output = apoE4_firstLastLaps(headscan_struct)
% Compare the first and last lap headscan metrics
    groups = {'apoE4', 'apoE3'};
    saveDir = '../figures/new headscans';
    FPS = 30;
    V_MAX = 6;

    for d = 1:4
        day_headscans = headscan_struct(strcmp([headscan_struct(:).day], num2str(d)));
        for s = 1:3
            session_headscans = day_headscans(strcmp([day_headscans(:).session], num2str(s)));
            figure(1); clf; hold on;

            % Make a matrix to store all of the values with 
            % dimensions: rats x group x first/last lap x metric number
            metric_data = nan(6, 2, 2, 6);
            for apo = 1:2
                group_headscans = session_headscans(strcmp({session_headscans(:).group}, groups(apo)));

                for h = 1:length(group_headscans)
                    session = group_headscans(h).headscan;

                    if(session.num_laps > 2)
                        % Use min and max to find the first and last laps
                        % (always from the end point)
                        firstLapEndTime = find(session.laps > min(session.laps)+1, 1);
                        lastLapStartTime = find(session.laps > max(session.laps)-1, 1);
                        
                        % First and last lap times
                        firstLapTime = (firstLapEndTime-1)/FPS;
                        lastLapTime = (length(session.laps)-lastLapStartTime)/FPS;
                        metric_data(h, apo, :, 1) = [firstLapTime, lastLapTime];
    
                        % Calculate the pause time
                        firstLapPause = sum(abs(session.dadt(1:firstLapEndTime)) < V_MAX)/FPS;
                        lastLapPause = sum(abs(session.dadt(lastLapStartTime:end)) < V_MAX)/FPS;
                        metric_data(h, apo, :, 2) = [firstLapPause, lastLapPause];
    
                        % Number of headscans
                        firstLapNumHeadscans = sum(session.filt_pss(:,1) < firstLapEndTime);
                        lastLapNumHeadscans = sum(session.filt_pss(:,1) > lastLapStartTime);
                        metric_data(h, apo, :, 3) = [firstLapNumHeadscans, lastLapNumHeadscans];
    
                        % Headscan Time
                        firstLapHeadscanTime = sum(session.filt_pss(session.filt_pss(:,1) < firstLapEndTime, 2) - session.filt_pss(session.filt_pss(:,1) < firstLapEndTime, 1))/FPS;
                        lastLapHeadscanTime = sum(session.filt_pss(session.filt_pss(:,1) > lastLapStartTime, 2) - session.filt_pss(session.filt_pss(:,1) > lastLapStartTime, 1))/FPS;
                        metric_data(h, apo, :, 4) = [firstLapHeadscanTime, lastLapHeadscanTime];
    
                        % Headscan ISI
                        firstLapHeadscanISI = firstLapTime / firstLapNumHeadscans;
                        lastLapHeadscanISI = lastLapTime / lastLapNumHeadscans;
                        metric_data(h, apo, :, 5) = [firstLapHeadscanISI, lastLapHeadscanISI];
    
                        % Scan Rate
                        firstLapScanRate = firstLapNumHeadscans / firstLapPause;
                        lastLapScanRate = lastLapNumHeadscans / lastLapPause;
                        metric_data(h, apo, :, 6) = [firstLapScanRate, lastLapScanRate];
                    end
                end
            end
            metric_mean = squeeze(mean(metric_data, 1, 'omitnan'));
            metric_std = squeeze(std(metric_data, [], 1, 'omitnan')) / sqrt(sum(~isnan(metric_data(:,apo,1,1))));

            for i = 1:6
                subplot(3,2,i); hold on;
                bar(metric_mean(:,:,i));

                % apoE4
                plot(0.85*ones(size(metric_data,1)), metric_data(:,1,1,i),'.k');
                plot(1.15*ones(size(metric_data,1)), metric_data(:,1,2,i),'.k');
                errorbar([0.85, 1.15], metric_mean(1,:,i), metric_std(1,:,i), 'vertical', 'Color','k');

                % apoE3
                plot(1.85*ones(size(metric_data,1)), metric_data(:,2,1,i),'.k');
                plot(2.15*ones(size(metric_data,1)), metric_data(:,2,2,i),'.k');
                errorbar([1.85, 2.15], metric_mean(2,:,i), metric_std(2,:,i), 'vertical', 'Color','k');

                legend({'First Lap','Last Lap'})
                xticks([1, 2])
                xticklabels({'apoE4','apoE3'})
            end
        
            title_pre = ['Day ' num2str(d) ' Session ' num2str(s)];
            subplot(3,2,1); title([title_pre ' Lap Time']); ylabel('Time [s]')
            subplot(3,2,2); title([title_pre ' Pause Time']); ylabel('Time [s]')
            subplot(3,2,3); title([title_pre ' Number of Headscans']); ylabel('Headscans')
            subplot(3,2,4); title([title_pre ' Headscan Time']); ylabel('Time [s]')
            subplot(3,2,5); title([title_pre ' Headscan ISI']); ylabel('ISI (Lap Time / # Headscans)')
            subplot(3,2,6); title([title_pre ' Scan Rate']); ylabel('# Headscans / Lap Pause Time')

            saveas(gcf, [saveDir filesep '20230703_Day' num2str(d) 'Session' num2str(s)], 'png')
        end
    end
end