function group = adCircTrack_nick_spikeTiming(group)
%ADCIRCTRACK_NICK_SPIKETRACK Summary of this function goes here
%   Detailed explanation goes here
    thetaBins = 250; % (Fs / FTheta) / 5 = (2000 / 8) / 5
    gammaBins = 45; % (Fs / FGamma) / 5
    figure(1); clf;

    saveDir = 'C:\Users\nrive\Projects\Colgin Lab\apoE4\figures\lfp_spikes';
    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                tetNums = group(g).rat(r).day(d).tetNums;
                tt = group(g).rat(r).day(d).thetaTet;
                cscFn = ['CSC' num2str(tetNums(tt)) '.ncs'];

                % Keep track of spike timings for these frequency ranges
                spikeTimeBins_theta = zeros(thetaBins, 1);
                spikeTimeBins_gamma = zeros(gammaBins, 1);
                spikeThetaGamma = zeros(thetaBins, gammaBins);

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
                    
                    counter = 0;
                    for i = 1:length(group(g).rat(r).day(d).begin(b).eegInds)
                        eegStart = group(g).rat(r).day(d).begin(b).eegInds(i,1);
                        eegStop = group(g).rat(r).day(d).begin(b).eegInds(i,2);

                        startTime = ts(eegStart);
                        stopTime = ts(eegStop);          

                        % Run cut2theta to get the timings to align spikes
                        % to
                        eegSig = bpLFP(eegStart:eegStop);
                        [time_theta,W,kInd] = cut2theta_nick(eegSig, lfpStruct.Fs);
                        [time_gamma,~,kInd_gamma] = cut2freq_nick(eegSig, 45, lfpStruct.Fs);
                            
                        % Add spikes to the times given
                        for u = 1:length(unit)
                            spkTms = unit(u).spkTms([unit(u).spkTms] > startTime & [unit(u).spkTms] < stopTime);
                            spkTms = floor((spkTms - startTime)*lfpStruct.Fs);

                            counter = counter+length(spkTms);

                            % Check if spkTms is empty first
                            if(~isempty(spkTms) && ~isempty(kInd))
                                % Check if it fits each window
                                for k = 1:size(kInd,1)
                                    kDiff_theta = floor(spkTms-kInd(k,1));
                                    kDiff = kDiff_theta(kDiff_theta > 0 & kDiff_theta <= length(time_theta));

                                    for s = 1:size(kDiff,1)
                                        thetaPhase = mod(kDiff(s),thetaBins)+1;
                                        gammaPhase = mod(kDiff(s),gammaBins)+1;
                                        spikeTimeBins_theta(thetaPhase) = spikeTimeBins_theta(thetaPhase)+1;
                                        spikeThetaGamma(thetaPhase, gammaPhase) = spikeThetaGamma(thetaPhase, gammaPhase) + 1;
                                    end
                                end

                                for k = 1:size(kInd_gamma,1)
                                    kDiff_gamma = floor(spkTms-kInd_gamma(k,1));
                                    kDiff = kDiff_gamma(kDiff_gamma > 0 & kDiff_gamma <= length(time_gamma));
                                    
                                    for s = 1:size(kDiff,1)
                                        gammaPhase = mod(kDiff(s),gammaBins)+1;
                                        spikeTimeBins_gamma(gammaPhase) = spikeTimeBins_gamma(gammaPhase)+1;
                                    end
                                end
                            end
                        end
                        
                    end
                end

                % Create heatmap and spike plot together
                figure(1); subplot(2,2,g); hold on;
                time_theta = linspace(0, 720, 2*thetaBins); % 1.44 = 360 (deg) / (2000Hz (Fs) / 8Hz (theta))
                smoothedSpikeTimeBins = movmean([spikeTimeBins_theta; spikeTimeBins_theta],10);
                normSSTB = smoothedSpikeTimeBins./sum(smoothedSpikeTimeBins);
                plot(time_theta, normSSTB)
                xlabel('Theta Phase (8Hz)')
                ylabel('Normalized Spike Probability (During running)')

                subplot(2,2,g+2); hold on;
                time_gamma = linspace(1, 360*2, 2*gammaBins);
                smoothedSpikeTimeBins = movmean([spikeTimeBins_gamma; spikeTimeBins_gamma],5);
                normSSTB = smoothedSpikeTimeBins./sum(smoothedSpikeTimeBins);
                plot(time_gamma, normSSTB)
                xlabel('Gamma Phase (45 Hz)')
                ylabel('Normalized Spike Probability (During running)')

                figure(); colormap(hot); imagesc(linspace(0, 360, thetaBins), linspace(-90, 270, gammaBins), spikeThetaGamma); colorbar 
                xlabel('Theta Phase'); ylabel('Gamma Phase')
                %subplot(2,1,2); colormap(hot); imagesc(thetaTime, 1:50, pow2db(Pxx));
                %saveas(gcf, [saveDir filesep '20230524_' group(g).name 'Day' num2str(d) '_thetaTet.png'])
            end
        end
    end

    figure(1);
    subplot(2,2,1); title('apoE4 Theta Modulation'); legend({'Day 1','Day 2','Day 3','Day 4','Day 5'})
    subplot(2,2,2); title('apoE3 Gamma Modulation'); legend({'Day 1','Day 2','Day 3'})
    subplot(2,2,3); title('apoE4 Theta Modulation'); legend({'Day 1','Day 2','Day 3','Day 4','Day 5'})
    subplot(2,2,4); title('apoE3 Gamma Modulation'); legend({'Day 1','Day 2','Day 3'})
end
