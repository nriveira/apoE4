function group = adCircTrack_nick_thetaTrace(group)
    saveDir = 'C:\Users\nrive\Projects\Colgin Lab\apoE4\figures\lfp_power';
    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                figure(1); clf; figure(2); clf;
                for b = 1:length(group(g).rat(r).day(d).begin)
                    figure(1); subplot(2,2,b); hold on;
                    title(['Begin ' num2str(b)])
                    xlabel('Normalized Time [s]')
                    ylabel('Z-scored Frequency')

                    samples = [];
                    thetaTime = group(g).rat(r).day(d).begin(b).thetaTime;

                    for c = 1:length(group(g).rat(r).day(d).begin(b).lapCuts)
                        thetaTrace = group(g).rat(r).day(d).begin(b).lapCuts(c).thetaCuts;
                        samples = [samples mean(thetaTrace,2)];
                        plot(thetaTime, mean(thetaTrace,2))
                    end

                    % Find xcorr values
                    corrVals = zeros(size(samples,2));
                    figure(2); subplot(2,2,b); hold on;
                    for i = 1:size(samples,2)
                        for j = 1:size(samples,2)
                            [xr, lags] = xcorr(samples(:,i), samples(:,j));
                            corrVals(i,j) = corr(samples(:,i), samples(:,j));
                            plot(lags, xr);
                        end
                    end
                    mean(corrVals, "all")
                    % Also include spectrogram
                    Pxx = get_wavelet_power(mean(samples,2), 2000, [1, 25], 6);
                    figure(3); colormap(hot); imagesc(thetaTime, 1:50, pow2db(Pxx)); colorbar; title([group(g).name ' Day ' num2str(d) ' Begin ' num2str(b)])
                    ylabel('Frequency')
                    xlabel('Normalized time [s]')
                    %saveas(gcf, [saveDir filesep '20230522_specTraceLap_' group(g).name 'Day' num2str(d) 'Begin' num2str(b) '.png'])
                end
                %figure(1); saveas(gcf, [saveDir filesep '20230522_rawTraceLap_' group(g).name 'Day' num2str(d) 'Begin' num2str(b) '.png'])
            end
        end
    end
end