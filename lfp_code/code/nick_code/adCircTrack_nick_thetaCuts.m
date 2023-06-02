function group = adCircTrack_nick_thetaCuts(group)
    %saveDir = 'C:\Users\nrive\Projects\Colgin Lab\apoE4\figures\lfp_power'; 
    for g = 1:length(group)
        samples = [];
        subplot(2,1,g); hold on;

        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                for b = 1:length(group(g).rat(r).day(d).begin)
                    thetaTime = group(g).rat(r).day(d).begin(b).thetaTime;
                    thetaCuts = group(g).rat(r).day(d).begin(b).thetaCuts;
                    samples = [samples thetaCuts];
                    plot(thetaTime, mean(thetaCuts,2))
                end
            end
        end

        % Waveform figure
        figure(1); title([group(g).name ' Average Waveform (per begin)'])
        ylim([-1.5 1.5])
        ylabel('Z-scored Amplitude (a.u.)')
        xlabel('Normalized time [s]')
        %saveas(gcf, [saveDir filesep '20230523_avgWaveformPerBegin_' group(g).name '.png'])

        Pxx = get_wavelet_power(mean(samples,2), 2000, [1, 25], 6);
        figure; colormap(hot); imagesc(thetaTime, 1:50, pow2db(Pxx)); colorbar; title([group(g).name ' Average Spectrogram (per begin)'])
        ylabel('Frequency')
        xlabel('Normalized time [s]')
        %saveas(gcf, [saveDir filesep '20230523_avgSpecPerBegin_' group(g).name '.png'])
    end
end