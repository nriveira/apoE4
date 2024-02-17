function group = adCircTrack_nick_spikeTiming(group)
%ADCIRCTRACK_NICK_SPIKETRACK Summary of this function goes here
    saveDir = 'C:\Users\nrive\Projects\Colgin Lab\apoE4\figures\updated spikes';

    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            thetaPhase_spkTms_Rat = zeros(1, 360);
            thetaPhase_cellNum_Rat = [];
            for d = 1:min(length(group(g).rat(r).day), 4)
                tetNums = group(g).rat(r).day(d).tetNums;
                tt = group(g).rat(r).day(d).thetaTet;
                cscFn = ['CSC' num2str(tt) '.ncs'];
                thetaPhase_spkTms_Day = zeros(1, 360);
                thetaPhase_cellNum_Day = [];
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
                        
                        for u = 1:length(unit)
                            spkTms = unit(u).spkTms(([unit(u).spkTms] > startTime) & ([unit(u).spkTms] < stopTime));
                            thetaPhase_cellNum_Day = [thetaPhase_cellNum_Day, (d*100)+u];
                            thetaPhase_cellNum_Rat = [thetaPhase_cellNum_Rat (d*100)+u];
                            for spk = 1:length(spkTms)
                                indx = ceil(theta_phase(find(spkTms(spk) < ts, 1)));
                                thetaPhase_spkTms_Day(indx) = thetaPhase_spkTms_Day(indx)+1;
                                thetaPhase_spkTms_Rat(indx) = thetaPhase_spkTms_Rat(indx)+1;
                            end
                        end
                    end
                end

                figure(1); clf; hold on;
                if(g == 1)
                    color = '#FF8800';
                else
                    color = '#48D1CC';
                end
                    
                spikes720 = zeros(90, 1);
                for s = 1:length(spikes720)/2
                    spikes720([s, s+45]) = sum(thetaPhase_spkTms_Day(8*(s-1)+1:min(8*s, length(thetaPhase_spkTms_Day))));
                end
                bar(1:8:720, spikes720, 'FaceColor', color)
                smoothed = conv(spikes720, ones(3, 1)./3, 'same');
                plot(9:8:712, smoothed(2:end-1), 'k', 'LineWidth', 2);
                ylabel('Number of Spikes')
                xlabel('Theta Phase')
                title(['Spikes per Theta Phase (' group(g).name ') Day ' num2str(d) ' (n=' num2str(length(unique(thetaPhase_cellNum_Day))) ')'])
                saveas(gcf, [saveDir filesep '20230717_SpikePhasePerDay' num2str(d) group(g).name], 'epsc')
                saveas(gcf, [saveDir filesep '20230717_SpikePhasePerDay' num2str(d) group(g).name], 'png')
            end

            figure(2); subplot(2,1,g); hold on;
            if(g == 1)
                color = '#FF8800';
            else
                color = '#48D1CC';
            end
            
            spikes720 = zeros(90, 1);
            for s = 1:length(spikes720)/2
                spikes720([s, s+45]) = sum(thetaPhase_spkTms_Rat(8*(s-1)+1:min(8*s, length(thetaPhase_spkTms_Rat))));
            end
            bar(1:8:720, spikes720, 'FaceColor', color)
            smoothed = conv(spikes720, ones(3, 1)./3, 'same');
            plot(9:8:712, smoothed(2:end-1), 'k', 'LineWidth', 2);
            ylabel('Number of Spikes')
            xlabel('Theta Phase')
            title([group(g).name ' n=' num2str(length(unique(thetaPhase_cellNum_Rat)))])
            saveas(gcf, [saveDir filesep '20230717_SpikePhasePerRat'], 'png')
            saveas(gcf, [saveDir filesep '20230717_SpikePhasePerRat'], 'epsc')
        end
    end
end
