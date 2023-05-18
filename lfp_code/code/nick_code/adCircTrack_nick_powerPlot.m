function group = adCircTrack_nick_powerPlot(group)
    colors = ['r','g','b','k'];
    for g = 1:length(group)
        figure(g); clf;
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                subplot(length(group(g).rat(r).day), 1, d); hold on;
                title(['Day ', num2str(d)])
                for b = 1:length(group(g).rat(r).day(d).begin)
                    % Get bootstrap confidence intervals of per begin power
                    % spectra
                    samples = pow2db(group(g).rat(r).day(d).begin(b).powerBegin');
                    power_freq = group(g).rat(r).day(d).begin(b).power_freq;
                    [ci, ~] = bootci(2000, {@mean, samples});

                    % After getting the bootstrap samples, plot them
                    lower = ci(1,:);
                    upper = ci(2,:);
                                        
                    patch([power_freq' fliplr(power_freq')], [lower fliplr(upper)], colors(b), 'EdgeColor','none', 'FaceAlpha', 0.3);
                    plot(power_freq, (upper+lower)/2, colors(b))
                    axis square; xlim([0, 25]); xline(4, '--');
                    xlabel('Frequency [Hz]')
                    ylabel('Power [dB]')
                end
                legend({'B1','', '','B2','','','B3','','','B4','',})
            end
        end
    end
end