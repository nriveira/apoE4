function output = nick_monacoAnalysis(group)
    PWELCH_WINDOW = 2048;
    PWELCH_PSD = 1025;
    FRAME_THRESH = 31;

    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                allDaysRunning = [];
                allDaysScanning = [];
                allDaysPausing = [];
                    
                for b = 1:length(group(g).rat(r).day(d).begin)
                    coords = group(g).rat(r).day(d).begin(b).coords;
                    timestamp = coords(:,1);
                    [dadt, filt_pss] = adCircTrack_x_findHeadscans(coords);

                    % Find running events (using dadt)
                    running_dadt = dadt > 6;
                    smoothed_dadt = conv(running_dadt, ones(floor(FRAME_THRESH/2),1), 'same');
                    running_event_index = find(smoothed_dadt > 1);
                    rii = 1;
                    prev = running_event_index(1);
                    running_events = [prev, 0];

                    for i = 2:length(running_event_index)
                        if(running_event_index(i) - prev == 1)
                            running_events(rii, 2) = running_event_index(i);
                        else
                            rii=rii+1;
                            running_events(rii, 1) = running_event_index(i);
                            running_events(rii, 2) = running_event_index(i);
                        end
                        prev = running_event_index(i);
                    end

                    %Filter out events less than 1 seconds
                    running_events(running_events(:,2)-running_events(:,1) < FRAME_THRESH,:) = [];
                    for re = 1:size(running_events,1)
                        running_events(re,:) = timestamp(running_events(re,:));
                    end

                    % Find Scanning events (using pss)
                    scanning_events = filt_pss;
                    scanning_events(scanning_events(:,2)-scanning_events(:,1) < FRAME_THRESH,:) = [];
                    for se = 1:size(scanning_events,1)
                        scanning_events(se, :) = timestamp(scanning_events(se, :));
                    end

                    % Find Paused events (Not moving fast and not in pss)
                    paused_event_index = find(dadt < 6);
                    paused_notPss = zeros(size(paused_event_index));
                    for i = 1:length(paused_event_index)
                        first_pss = find(filt_pss(:,1) > paused_event_index(i), 1);
                        last_pss = find(filt_pss(:,2) < paused_event_index(i), 1);
                        if(first_pss ~= last_pss) % If this is true it must be in the pss, so dont include
                            paused_notPss(i) = 1;
                        end
                    end
                    
                    pii = 1;
                    prev = paused_event_index(1);
                    paused_events = [prev, 0];

                    for i = 2:length(paused_event_index)
                        if(paused_event_index(i) - prev == 1)
                            paused_events(pii, 2) = paused_event_index(i);
                        else
                            pii=pii+1;
                            paused_events(pii, 1) = paused_event_index(i);
                            paused_events(pii, 2) = paused_event_index(i);
                        end
                        prev = paused_event_index(i);
                    end
                    paused_events(paused_events(:,2)-paused_events(:,1) < FRAME_THRESH,:) = [];
                    for pe = 1:size(paused_events,1)
                        paused_events(pe, :) = timestamp(paused_events(pe, :));
                    end

                    % Load in the LFP data to get the PSD
                    theta_tet = group(g).rat(r).day(d).thetaTet;
                    lfp_dir = group(g).rat(r).day(d).begin(b).dir;
                    cd(lfp_dir)

                    cscFn = ['CSC' num2str(theta_tet) '.ncs'];
                    lfpStruct = read_in_lfp(cscFn);

                    % Bandpass filter lfp data
                    bpLFP = bandpass(zscore(lfpStruct.data), [1 475], lfpStruct.Fs);
                    [pxx_norm, f] = pwelch(bpLFP, PWELCH_WINDOW, [], [], lfpStruct.Fs);
                    lfp_time = lfpStruct.ts;

                    % Align all of the running times with the lfp indices
                    for ii = 1:size(running_events,1)
                        startTimeRunning = running_events(ii,1);
                        stopTimeRunning = running_events(ii,2);
                        [~,startInd] = min(abs(lfp_time-startTimeRunning));
                        [~,stopInd] = min(abs(lfp_time-stopTimeRunning));
                        running_events(ii,:) = [startInd, stopInd];
                    end
                    running_lfp_collated = zeros(size(running_events,1), PWELCH_PSD);
                    for ii = 1:size(running_events,1)
                        startInd = fix(running_events(ii,1));
                        stopInd = fix(running_events(ii,2));
                        lfp_running = bpLFP(startInd:stopInd);
                        [pxx_running, f_running] = pwelch(lfp_running, PWELCH_WINDOW, [], [], lfpStruct.Fs);
                        running_lfp_collated(ii, :) = pxx_running;
                    end
                    group(g).rat(r).day(d).begin(b).running_PSD = running_lfp_collated;
                    allDaysRunning = [allDaysRunning; running_lfp_collated];
                    
                    % Align all of the scanning times with the lfp indices
                    for ii = 1:size(scanning_events,1)
                        startTimeRunning = scanning_events(ii,:);
                        [~,startInd] = min(abs(lfp_time-startTimeRunning(1)));
                        [~,stopInd] = min(abs(lfp_time-startTimeRunning(2)));
                        scanning_events(ii,:) = [startInd, stopInd];
                    end
                    scanning_lfp_collated = zeros(size(scanning_events,1), PWELCH_PSD);
                    for ii = 1:size(scanning_events,1)
                        startInd = fix(scanning_events(ii,1));
                        stopInd = fix(scanning_events(ii,2));
                        lfp_scanning = bpLFP(startInd:stopInd);
                        [pxx_scanning, f_scanning] = pwelch(lfp_scanning, PWELCH_WINDOW, [], [], lfpStruct.Fs);
                        scanning_lfp_collated(ii, :) = (pxx_scanning)';
                    end
                    group(g).rat(r).day(d).begin(b).scanning_PSD = scanning_lfp_collated;
                    allDaysScanning = [allDaysScanning; scanning_lfp_collated];

                    % Align all of the pause times with the lfp indices
                    for ii = 1:size(paused_events,1)
                        startTimeRunning = paused_events(ii,:);
                        [~,startInd] = min(abs(lfp_time-startTimeRunning(1)));
                        [~,stopInd] = min(abs(lfp_time-startTimeRunning(2)));
                        paused_events(ii,:) = [startInd, stopInd];
                    end
                    paused_lfp_collated = zeros(size(paused_events,1), PWELCH_PSD);
                    for ii = 1:size(paused_events,1)
                        startInd = fix(paused_events(ii,1));
                        stopInd = fix(paused_events(ii,2));
                        lfp_paused = bpLFP(startInd:stopInd);
                        [pxx_paused, f_paused] = pwelch(lfp_paused, PWELCH_WINDOW, [], [], lfpStruct.Fs);
                        paused_lfp_collated(ii, :) = (pxx_paused)';
                    end
                    group(g).rat(r).day(d).begin(b).paused_PSD = paused_lfp_collated;
                    allDaysPausing = [allDaysPausing; paused_lfp_collated];
                end

                allDaysRunning = (allDaysRunning./std(allDaysRunning, [], 1));
                allDaysScanning = (allDaysScanning./std(allDaysScanning, [], 1));
                allDaysPausing = (allDaysPausing./std(allDaysPausing, [], 1));

                % Plot all from a day
                figure(g); hold on;
                plot(f_running, mean(allDaysRunning-mean(allDaysRunning), 1), 'k')
                plot(f_scanning, mean(allDaysScanning-mean(allDaysScanning), 1), 'b')
                plot(f_paused, mean(allDaysPausing-mean(allDaysPausing), 1), 'g')
%                 plot_conf(allDaysRunning, f_running, 'k',['Running Day ' num2str(d)])
%                 plot_conf(allDaysScanning, f_running, 'b',['Scanning Day ' num2str(d)])
%                 plot_conf(allDaysPausing, f_running, 'g',['Paused Day ' num2str(d)])
                xlim([4, 12])
                legend()
            end
        end
    end
end

function output = plot_conf(eventMatrix, f, color, displayName)
    event_95conf = std(eventMatrix, [], 1);
    event_mean = mean(eventMatrix, 1);
    y1 = event_mean + event_95conf;
    y2 = event_mean - event_95conf;
    patch([f', fliplr(f')], [y1, fliplr(y2)], color, 'DisplayName', displayName, 'FaceAlpha', 0.5, 'EdgeAlpha', 0)
end