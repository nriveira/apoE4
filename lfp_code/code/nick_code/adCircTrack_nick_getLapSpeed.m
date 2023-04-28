% Add radial velocity (in cm) to the group struct
function group = adCircTrack_nick_getLapSpeed(group)
conversion_factor = (pi * 50) / (180*30); % Converting to cm/s given track radius 

for g = 1:2
%for g = 1 %this is only for apoE4 rat
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            tic
            for b = 1:length(group(g).rat(r).day(d).begin)
                fprintf('\t\t\tBegin %d\n', b);
                velocity = diff(group(g).rat(r).day(d).begin(b).radPos(:,2))./conversion_factor;
                velocity(abs(velocity) > 90) = 0;
                timestamps = group(g).rat(r).day(d).begin(b).radPos(2:end,1);
                fs = 1/(timestamps(2)-timestamps(1));

                winLength = floor(fs/2);
                velocity = abs(smoothdata(velocity, "gaussian", winLength));

                velocityBin = zeros(floor(length(velocity)/winLength),2);
                % Aggregate velocities into 0.5 second non-overlapping time bins
                for i = 1:length(velocityBin)
                    indices = winLength*(i-1)+1:winLength*i;
                    velocityBin(i,1) = mean(timestamps(indices));
                    velocityBin(i,2) = mean(velocity(indices));
                end

                group(g).rat(r).day(d).begin(b).velocityBin = velocityBin;
                
                % Access the EEG data to get power values
                cscFn = ['CSC' num2str(group(g).rat(r).day(d).thetaTet) '.ncs'];
                fprintf('\t\t\t\tFiltering LFP for Tetrode #%d\n', group(g).rat(r).day(d).thetaTet);
                lfpStruct = read_in_lfp([group(g).rat(r).day(d).begin(b).dir filesep cscFn]);
                eeg_Z = bandpass(lfpStruct.data, [1 250], lfpStruct.Fs);
                
                % For each timestep found in the velocity bin, calculate
                % the slow gamma and fast gamma powers
                gamma_wavelet = zeros(size(velocityBin));

                for i = 1:length(gamma_wavelet)
                    start_ind = find(lfpStruct.ts > velocityBin(i,1), 1);
                    eeg_window = max((-lfpStruct.Fs/4)+start_ind, 1):min((lfpStruct.Fs/4)+start_ind, length(eeg_Z));

                    eeg = eeg_Z(eeg_window);
                    slow_g = sum(get_wavelet_power(eeg, lfpStruct.Fs, [25 55], 6), 'all');
                    fast_g = sum(get_wavelet_power(eeg, lfpStruct.Fs, [70 130], 6), 'all');
                    gamma_wavelet(i, :) = [slow_g, fast_g];
                end

                gamma_wavelet_zscore = zscore(gamma_wavelet);
                group(g).rat(r).day(d).begin(b).gamma_wavelet = gamma_wavelet_zscore;

                figure(1); hold on;
                plot(log2(velocityBin(:,2)), gamma_wavelet_zscore(:,1), '.');

                figure(2); hold on;
                plot(log2(velocityBin(:,2)), gamma_wavelet_zscore(:,2), '.');
            end

            figure(1); xlim([-2 6]); ylim([-3 6]); title([group(g).name ' Day ' num2str(d) ' Slow Gamma'])
            xlabel('log2(Speed) (cm/s)'); ylabel('Slow Gamma Power [stdev from mean]');
            save_name = ['.\figures\20230414_' group(g).name 'd' num2str(d) '_slowgamma.png'];
            saveas(gcf, save_name); clf;
            figure(2); xlim([-2 6]); ylim([-3 6]); title([group(g).name ' Day ' num2str(d) ' Fast Gamma'])
            xlabel('log2(Speed) (cm/s)'); ylabel('Fast Gamma Power [stdev from mean]');
            save_name = ['.\figures\20230414_' group(g).name 'd' num2str(d) '_fastgamma.png'];
            saveas(gcf, save_name); clf;
            toc
        end
    end
end

end