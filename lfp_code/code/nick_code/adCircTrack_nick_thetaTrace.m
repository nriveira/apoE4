function group = adCircTrack_nick_thetaTrace(group)
    saveDir = 'C:\Users\nrive\Projects\Colgin Lab\apoE4\figures\lfp_power';
    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                figure(1); clf;
                for b = 1:length(group(g).rat(r).day(d).begin)
                    figure(2); clf;
                    thetaTime = group(g).rat(r).day(d).begin(b).thetaTime;
                    for c = 1:length(group(g).rat(r).day(d).begin(b).lapCuts)
                        thetaTrace = group(g).rat(r).day(d).begin(b).lapCuts(c).thetaCuts;
                        figure(1); subplot(2,2,b); hold on; plot(thetaTime, mean(thetaTrace,2))
                        title([group(g).name ' Waveform per Lap Begin ' num2str(b)])
                        xlabel('Normalized Time [s]')
                        ylabel('Z-scored Frequency')

                        % Also get spec per lap
                        Pxx = get_wavelet_power(mean(thetaTrace,2), 2000, [1, 25], 6);
                        figure(2); colormap(hot); imagesc(thetaTime, 1:50, pow2db(Pxx)); colorbar; title([group(g).name ' d' num2str(d) 'b' num2str(b) ' lap' num2str(c)])
                        ylabel('Frequency')
                        xlabel('Normalized time [s]')
                        figure(2); saveas(gcf, [saveDir filesep '20230523_specTraceLap_' group(g).name 'Day' num2str(d) 'Begin' num2str(b) 'Lap' num2str(c) '.png'])
                    end
                end
                figure(1); saveas(gcf, [saveDir filesep '20230523_avgTraceLap_' group(g).name 'Day' num2str(d) '.png'])
            end
        end
    end
end