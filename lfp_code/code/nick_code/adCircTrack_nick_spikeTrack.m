function group = adCircTrack_nick_spikeTrack(group)
%ADCIRCTRACK_NICK_SPIKETRACK Summary of this function goes here
%   Detailed explanation goes here

    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                tetNums = group(g).rat(r).day(d).tetNums;
                tt = group(g).rat(r).day(d).thetaTet;
                cscFn = ['CSC' num2str(tetNums(tt)) '.ncs'];
                spikeTimeBins = zeros(1201, 1);

                for b = 1:length(group(g).rat(r).day(d).begin)
                    cd(group(g).rat(r).day(d).begin(b).dir)
                    % Load in lfp data and eeg indices
                    lfpStruct = read_in_lfp(cscFn);
                    
                    ts = lfpStruct.ts;
                    bpLFP = bandpass(zscore(lfpStruct.data), [0.5 200], lfpStruct.Fs);
                    unit = group(g).rat(r).day(d).begin(b).unit;

                    unitInd = reshape([unit.ID],2,[]);
                    unitInd(2,:) = [];
                    unit = unit(unitInd == tetNums(tt));
                    

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
                            spkTms = unit(u).spkTms;%([unit(u).spkTms] > startTime & [unit(u).spkTms] < stopTime);
                            spkTms = floor((spkTms - startTime)*lfpStruct.Fs);

                            % Check if spkTms is empty first
                            if(~isempty(spkTms))
                                % Go through each spike
                                for s = 1:length(spkTms)
                                    % Check if it fits each window
                                    for k = 1:length(kInd)
                                        kDiff = floor(kInd(k,1)-spkTms(s));
                                        if(kDiff > 0 & kDiff <= 1200)
                                            spikeTimeBins(kDiff) = spikeTimeBins(kDiff)+1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                figure(); 
                bar(time, movmean(spikeTimeBins,7),5)
                xlabel('Normalized Time [s]')
                ylabel('Number of spikes 7-frame moving average (all tets)')
                title([group(g).name ' Day ' num2str(d) ' (n = ' num2str(sum(spikeTimeBins)) ')'])
                %group(g).rat(r).day(d).spikeTimeBins = spikeTimeBins;
            end
        end
    end
end

