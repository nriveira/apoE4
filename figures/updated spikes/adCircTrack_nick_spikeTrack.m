function group = adCircTrack_nick_spikeTrack(group)
%ADCIRCTRACK_NICK_SPIKETRACK Summary of this function goes here
%   Detailed explanation goes here
    colors = {'#48D1CC', '#FF8800'};
    saveDir = 'C:\Users\nrive\Projects\Colgin Lab\apoE4\figures\updated spikes';
    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            spikePhase = zeros(360,1);
            for d = 1:min(length(group(g).rat(r).day), 4)
                tt = group(g).rat(r).day(d).thetaTet;
                cscFn = ['CSC' num2str(tt) '.ncs'];
                for b = 1:length(group(g).rat(r).day(d).begin)
                    cd(group(g).rat(r).day(d).begin(b).dir)
                    % Load in lfp data and eeg indices
                    lfpStruct = read_in_lfp(cscFn);
                    
                    ts = lfpStruct.ts;
                    bpLFP = bandpass(zscore(lfpStruct.data), [1 250], lfpStruct.Fs);
                    unit = group(g).rat(r).day(d).begin(b).unit;

                    thetaLFP = bandpass(zscore(lfpStruct.data), [6 12], lfpStruct.Fs);
                    theta_phase = rad2deg(angle(hilbert(thetaLFP)))+180;

%                     unitInd = reshape([unit.ID],2,[]);
%                     unitInd(2,:) = [];
%                     unit = unit(unitInd == tt);
                    
                    radPos = group(g).rat(r).day(d).begin(b).radPos;
                    coords = group(g).rat(r).day(d).begin(b).coords;
                    velocity = get_runspeed(coords);
                    velocity = smooth_runspeed(velocity);

                    % Find running events (> 5cm/s) If running events are
                    % within 0.5 seconds of each other and the rat is
                    % running for more than half of the time
                    running = velocity(:,2) > 5;
                    
                    is_running = diff(running);
                    running_event = [];
                    running_vector = zeros(length(velocity),1);

                    % If the end flag comes first, then the first frame
                    % is a running frame
                    if(find(is_running == 1, 1) > find(is_running == -1,1))
                        running_vector(1) = 1;
                    end
                    running_vector(2:end) = is_running;
                    run_starts = find(running_vector == 1);
                    run_ends = find(running_vector == -1);
                    if(length(run_starts) ~= length(run_ends))
                        run_ends = [run_ends; length(velocity)];
                    end

                    % Eliminate all run events less than 1 second
                    run_threshold = (run_ends - run_starts) < 30;
                    run_starts(run_threshold) = [];
                    run_ends(run_threshold) = [];

                    for rs = 1:length(run_starts)
                        startTime = find(radPos(run_starts(rs),1) < ts, 1);
                        stopTime = find(radPos(run_ends(rs),1) < ts, 1);
                        eegSignal = bpLFP(startTime:stopTime);

                        [time, ~, kInd] = cut2theta_nick(eegSignal, lfpStruct.Fs);
                        
                        % For each cell
                        for u = 1:length(unit)
                            spkTms = unit(u).spkTms(([unit(u).spkTms] > startTime) & ([unit(u).spkTms] < stopTime));
                            % For each spike
                            for spk = 1:length(spkTms)
                                indx = find(spkTms(spk) < ts, 1);
                                % For each kInd
                                for k = 1:length(kInd)
                                    if(abs((indx-startTime)-kInd(k)) < 0.3*lfpStruct.Fs)
                                        spkPhase = ceil(theta_phase(indx));
                                        spikePhase(spkPhase) = spikePhase(spkPhase)+1;
                                    end
                                end
                            end
                        end
                    end
                end
            end

            spikePhase30 = zeros(12, 1);

            for i = 1:length(spikePhase30)
                spikePhase30(i) = sum(spikePhase(30*(i-1)+1:30*i));
            end 

            subplot(2,1,g); 
            b = bar(1:30:720, [spikePhase30; spikePhase30]);
            b.FaceColor = colors{g};
            title(['Spikes During Theta Mode (n=' num2str(sum(spikePhase)) ')'])
            xlabel('Theta Phase')
            ylabel('Number of Spikes')
        end
    end
    saveas(gcf, [saveDir filesep '20230717_spikesDuringTheta'], 'png')
    saveas(gcf, [saveDir filesep '20230717_spikesDuringTheta'], 'epsc')
end
