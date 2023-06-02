function group = adCircTrack_nick_power(group)
    % Tag laps first (Only run once)
    %group = adCircTrack_3_tagLaps_jayanth(group);
    FPS = 30;
    TRACK_RAD = 50;
    CM_THRESH = 5;
    RUN_THRESH = 1 * FPS;

    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                tetNums = group(g).rat(r).day(d).tetNums;
                tt = group(g).rat(r).day(d).thetaTet;
                cscFn = ['CSC' num2str(tetNums(tt)) '.ncs'];

                for b = 1:length(group(g).rat(r).day(d).begin)
                    cd(group(g).rat(r).day(d).begin(b).dir)
                    %fprintf('\t\t\t\tFiltering LFP for Tetrode #%d\n', tetNums(tt));

                    % Load in position and lfp data and lap indices
                    lfpStruct = read_in_lfp(cscFn);
                    [~,posy,posx] = LoadCircPos('VT1.nvt');
                    radPos = circpos(posx,posy);
                    lfpTime = group(g).rat(r).day(d).begin(b).radPos(:,1);

                    % Create Lap Numbers array
                    lapInds = group(g).rat(r).day(d).begin(b).lapInds;
                    lapNum = zeros(size(lfpTime));
                    totalLaps = length(lapInds);
                    for tl = 1:totalLaps
                        lapNum(lapInds(tl,1):lapInds(tl,2)) = tl;
                    end

                    % Bandpass filter lfp data
                    bpLFP = bandpass(zscore(lfpStruct.data), [0.5 200], lfpStruct.Fs);

                    % Clean up velocity measurements
                    degVelo = zeros(size(radPos));
                    degVelo(2:end) = TRACK_RAD*FPS*diff(deg2rad(radPos));
                    degVelo(abs(degVelo) > 100) = 0;
                    degVelo = abs(movmean(degVelo,7));

                    % Find all of the times moving over 6 cm/s and add one
                    % second buffer to each side
                    overThresh = degVelo > CM_THRESH; 
                    overThreshBuffer = conv(overThresh, ones(FPS/2, 1),'same') > 2;
%                     fprintf('Time in theta mode: %1.4f\n', sum(overThreshBuffer)/length(overThreshBuffer))
                    running = diff(overThreshBuffer);

                    %  For each instance of moving over 6 cm, run cut2theta
                    % if running for over 0.5 seconds
                    running_start = find(running == 1);
                    running_end = find(running == -1);

                    if(running_start(1) > running_end(1))
                        running_start = [1, running_start];
                    end

                    if(length(running_start) > length(running_end))
                        running_end = [running_end, length(running)];
                    end

                    running_length = running_end - running_start;
                    running_start(running_length < RUN_THRESH) = [];
                    running_end(running_length < RUN_THRESH) = [];
                    
                    fprintf('Seconds in theta mode: %3.3f\n', sum([running_end - running_start])/FPS)

                    lapTrack = lapNum(running_start);
                    runningInd = [running_start', running_end', lapTrack];
                    eegStarts = [];
                    eegStops = [];
                    thetaCutsBegin = [];
                    totalPowerBegin = [];
                    thetaPowerBegin = [];

                    % Loop through those times and find theta sequences
                    for tl = 1:totalLaps
                        curLap = runningInd(runningInd(:,3) == tl,:);
                        thetaCutsLap = [];
                        totalPowerLap = [];
                        thetaPowerLap = [];
                        for i = 1:size(curLap,1)
                            eegStart = find(lfpStruct.ts > lfpTime(curLap(i,1)), 1);
                            eegStop = find(lfpStruct.ts > lfpTime(curLap(i,2)), 1);

                            [time, W, ~] = cut2theta_nick(bpLFP(eegStart:eegStop), lfpStruct.Fs);
                            [pxx, f] = pwelch(bpLFP(eegStart:eegStop), hamming(1024), [], 1024, lfpStruct.Fs);
                            thetaPower = sum(pxx(f > 7 & f < 13));

                            % Also find the spikes that happened during
                            % that time to add to spike plot

                            % Get waveform measures
                            thetaCutsBegin = [thetaCutsBegin W];
                            thetaCutsLap = [thetaCutsLap W];

                            totalPowerBegin = [totalPowerBegin pxx];
                            totalPowerLap = [totalPowerLap pxx];

                            thetaPowerBegin = [thetaPowerBegin thetaPower];
                            thetaPowerLap = [thetaPowerLap thetaPower];

                            eegStarts = [eegStarts; eegStart];
                            eegStops = [eegStops; eegStop];

                        end
                        group(g).rat(r).day(d).begin(b).lapCuts(tl).thetaCuts = thetaCutsLap;
                        group(g).rat(r).day(d).begin(b).lapCuts(tl).totalPowerLap = totalPowerLap;
                        group(g).rat(r).day(d).begin(b).lapCuts(tl).thetaPowerLap = thetaPowerLap;
                    end
                    group(g).rat(r).day(d).begin(b).thetaCuts = thetaCutsBegin;
                    group(g).rat(r).day(d).begin(b).powerBegin = totalPowerBegin;
                    group(g).rat(r).day(d).begin(b).thetaPowerBegin = thetaPowerBegin;
                    group(g).rat(r).day(d).begin(b).eegInds = [eegStarts eegStops];
                    group(g).rat(r).day(d).begin(b).power_freq = f;
                    group(g).rat(r).day(d).begin(b).thetaTime = time;

                    %figure(g); subplot(5,1,d); hold on; plot(time, mean(thetaCutsBegin,2)); title(['Day' int2str(d) ' Number of Theta cuts: ' int2str(totalCuts)]);
                end
            end
        end
    end
end