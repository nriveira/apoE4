function group = adCircTrack_nick_spikeTrack(group)
%ADCIRCTRACK_NICK_SPIKETRACK Summary of this function goes here
%   Detailed explanation goes here
    colors = ['b', 'k'];
    saveDir = 'C:\Users\nrive\Projects\Colgin Lab\apoE4\figures\lfp_e3e4';
    for g = 1:length(group)
        subplot(2,1,g); hold on;
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                tetNums = group(g).rat(r).day(d).tetNums;
                tt = group(g).rat(r).day(d).thetaTet;
                cscFn = ['CSC' num2str(tt) '.ncs'];
                spikeTimeBins = zeros(250, 1);
                samples = [];

                for b = 1:length(group(g).rat(r).day(d).begin)
                    cd(group(g).rat(r).day(d).begin(b).dir)
                    % Load in lfp data and eeg indices
                    lfpStruct = read_in_lfp(cscFn);
                    
                    ts = lfpStruct.ts;
                    bpLFP = bandpass(zscore(lfpStruct.data), [0.5 200], lfpStruct.Fs);
                    unit = group(g).rat(r).day(d).begin(b).unit;

                    unitInd = reshape([unit.ID],2,[]);
                    unitInd(2,:) = [];
                    unit = unit(unitInd == tt);
                    

                    for i = 1:length(group(g).rat(r).day(d).begin(b).eegInds)
                        eegStart = group(g).rat(r).day(d).begin(b).eegInds(i,1);
                        eegStop = group(g).rat(r).day(d).begin(b).eegInds(i,2);

                        startTime = ts(eegStart);
                        stopTime = ts(eegStop);          

                        % Run cut2theta to get the timings to align spikes
                        % to
                        eegSig = bpLFP(eegStart:eegStop);
                        [time,W,kInd] = cut2theta_nick(eegSig, lfpStruct.Fs);
                            
                        % Add spikes to the times given
                        for u = 1:length(unit)
                            spkTms = unit(u).spkTms([unit(u).spkTms] > startTime & [unit(u).spkTms] < stopTime);
                            spkTms = floor((spkTms - startTime)*lfpStruct.Fs);

                            % Check if spkTms is empty first
                            if(~isempty(spkTms))
                                % Go through each spike
                                for s = 1:length(spkTms)
                                    % Check if it fits each window
                                    for k = 1:size(kInd,1)
                                        kDiff = floor(spkTms(s)-kInd(k,1));
                                        if(kDiff > 0 & kDiff <= 1200)
                                            spikeTimeBins(mod(kDiff,250)+1) = spikeTimeBins(mod(kDiff,250)+1)+1;
                                        end
                                    end
                                end
                            end
                        end
                    end

                    % Collect the theta traces to do the spec plot
                    thetaTime = group(g).rat(r).day(d).begin(b).thetaTime;
                    thetaCuts = group(g).rat(r).day(d).begin(b).thetaCuts;
                    samples = [samples thetaCuts];

                end

                % Create heatmap and spike plot together
                %fig_index = fig_index+1; 
                time = 1.44:1.44:720; % 1.44 = 360 (deg) / (2000Hz (Fs) / 8Hz (theta))
                smoothedSpikeTimeBins = movmean([spikeTimeBins; spikeTimeBins],10);
                normSSTB = smoothedSpikeTimeBins./sum(smoothedSpikeTimeBins);
                plot(time, normSSTB, 'DisplayName', [group(g).name ' Day ' num2str(d) ' (n = ' num2str(sum(spikeTimeBins)) ')'])
                xlabel('Theta Phase (8Hz)')
                ylabel('Posterior Spike Probability (During running)')
                title([group(g).name ' Phase coupling'])
                legend()
            end
        end
    end
end

