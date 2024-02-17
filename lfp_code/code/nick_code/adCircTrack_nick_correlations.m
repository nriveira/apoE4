function group = adCircTrack_nick_correlations(group)
%ADCIRCTRACK_NICK_SPIKETRACK Summary of this function goes here
    saveDir = 'C:\Users\nrive\Projects\Colgin Lab\apoE4\figures\updated spikes';

    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            for d = 1:min(length(group(g).rat(r).day), 4)
                tt = group(g).rat(r).day(d).thetaTet;
                cscFn = ['CSC' num2str(tt) '.ncs'];
                for b = 1:length(group(g).rat(r).day(d).begin)
                    cd(group(g).rat(r).day(d).begin(b).dir)

                    % Load in lfp eeg indices
                    lfpStruct = read_in_lfp(cscFn);
                    ts = lfpStruct.ts;                  
                    unit = group(g).rat(r).day(d).begin(b).unit;

%                     unitInd = reshape([unit.ID],2,[]);
%                     unitInd(2,:) = [];
%                     unit = unit(unitInd == tt);
                    new_size = floor(length(ts)/20);
                    spikes = zeros(new_size, length(unit));
                    for u = 1:length(unit)
                        for s = 1:length(unit(u).spkTms)
                            indx = ceil(find(unit(u).spkTms(s) < ts, 1)/20);
                            spikes(indx,u) = spikes(indx,u)+1;
                        end
                    end

                    figure(1); subplot(4,1,b); hold on;
                    title(['Begin ' num2str(b)])
                    for i = 1:length(unit)
                        [corr_r, lags] = xcorr(spikes(:,i), spikes(:,i), 30);
                        plot(0.01*lags, corr_r./sqrt(sum(spikes(:,i))), 'DisplayName', ['Unit ' num2str(i)]);
                    end
                end
            end
            clf;
        end
    end
end
