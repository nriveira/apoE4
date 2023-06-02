function group = adCircTrack_nick_thetaTrace(group)
    saveDir = 'C:\Users\nrive\Projects\Colgin Lab\apoE4\figures\lfp_power';
    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                figure(1); clf;
                for b = 1:length(group(g).rat(r).day(d).begin)

                    thetaTime = group(g).rat(r).day(d).begin(b).thetaTime;
                    for c = 1:length(group(g).rat(r).day(d).begin(b).lapCuts)
                        thetaTrace = group(g).rat(r).day(d).begin(b).lapCuts(c).thetaCuts;
                        figure(1); subplot(8,1,2*(b-1)+1); hold on; plot(thetaTime+(c), mean(thetaTrace,2), 'k')
                        title([group(g).name ' Waveform per Lap Begin ' num2str(b)])
                        xlabel('Normalized Time [s]')
                        ylabel('Z-scored Amplitude')

                        % Also get spec per lap
                        Pxx = get_wavelet_power(mean(thetaTrace,2), 2000, [1, 25], 6);
                        subplot(8,1,2*(b)); hold on; colormap(hot); imagesc(thetaTime+(c), 1:50, pow2db(Pxx)); colorbar; title([group(g).name ' Day ' num2str(d) ' Begin ' num2str(b)])
                        ylabel('Frequency')
                        xlabel('Lap # [-0.3, +0.3] window from onset')
                    end
                end
                figure(1); saveas(gcf, [saveDir filesep '20230529_waveSpec_' group(g).name 'Day' num2str(d) '.png'])
            end
        end
    end
end