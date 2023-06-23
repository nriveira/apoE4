function group = adCircTrack_nick_spikeLocating(group)
%ADCIRCTRACK_NICK_SPIKETRACK Summary of this function goes here
%   Detailed explanation goes here
    saveDir = 'C:\Users\nrive\Projects\Colgin Lab\apoE4\figures\lfp_e3e4';
    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                tetNums = group(g).rat(r).day(d).tetNums;
                tt = group(g).rat(r).day(d).thetaTet;
                cscFn = ['CSC' num2str(tt) '.ncs'];
                figure(d); hold on;

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
                    spikeStruct = [];

                    % For every lap, look at the spike rate as a function
                    % of theta phase
                    for lap = 1:length(group(g).rat(r).day(d).begin(b).lapTms)
                        toi = group(g).rat(r).day(d).begin(b).lapTms(lap,:);
                        eegInds = group(g).rat(r).day(d).begin(b).eegInds;
                        eegInds_toi = eegInds(eegInds(:,2) < find(lfpStruct.ts > toi(2), 1),:);
                        % Go through each running bout and only include
                        % ones that are in the time of interest
                        for i = 1:length(eegInds_toi)
                            eegStart = eegInds_toi(i,1);
                            eegStop = eegInds_toi(i,2);

                            startTime = ts(eegStart);
                            stopTime = ts(eegStop);

                            eegSig = bpLFP(eegStart:eegStop);
                            [time_theta, spikeStruct(lap).W, spikeStruct(lap).kInd] = cut2theta_nick(eegSig, lfpStruct.Fs);
                        end

                        for u = 1:length(unit)
                            spkTms = unit(u).spkTms([unit(u).spkTms] > startTime & [unit(u).spkTms] < stopTime);
                            spkTms = floor((spkTms - startTime)*lfpStruct.Fs);
                            kSpikes = [];
                            for k = 1:length(spikeStruct(lap).kInd)
                                kDiff = floor(spkTms - spikeStruct(lap).kInd(k,1));
                                kSpikes = [kSpikes; kDiff(kDiff > 0 & kDiff <  length(time_theta))];
                            end
                            spikeStruct(lap).spikes = kSpikes;
                        end
                    end
                    
                    if(g == 1)
                        subplot(4,2,2*(b-1)+1); hold on;
                    else
                        subplot(4,2,2*b); hold on;
                    end
                    plot_session(spikeStruct, time_theta);
                    title([group(g).name ' Day ' num2str(d) ' Begin ' num2str(b)])
                end %begin
            end %day
        end %rat
    end %group
end

function fig = plot_session(spikeStruct, time_theta)
    for i = 1:length(spikeStruct)
        if(~isempty(spikeStruct(i).W))
            plot(time_theta+i, mean(spikeStruct(i).W,2),'k')
            if(~isempty(spikeStruct(i).spikes))
                % Get location of each spike (relative to start index)
                spike_loc = histcounts(spikeStruct(i).spikes,1201) > 0;
                plot(time_theta(spike_loc)+i, mean(spikeStruct(i).W(spike_loc,:),2), 'b*')
            end
        end
    end
end