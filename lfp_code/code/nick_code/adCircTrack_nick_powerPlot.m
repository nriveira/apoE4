function group = adCircTrack_nick_powerPlot(group)
    % adCircTrack_nick_powerPlot
    % Plots the bootstrap estimated power across each begin for apoE3 and
    % apoE4 rats
    colors = ['b','k'];
    figure(1); clf;
    saveDir = 'C:\Users\nrive\Projects\Colgin Lab\apoE4\figures\lfp_e3e4'; 
    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            for d = 1:3%length(group(g).rat(r).day)
                % Create new figure per day and give title
                subplot(1,3,d); hold on; title(['apoE4 vs apoE3 bootstrap PSD during running Day ', num2str(d)])
                for b = 1:length(group(g).rat(r).day(d).begin)
                    % Get bootstrap confidence intervals of per begin power
                    % spectra
                    samples = pow2db(group(g).rat(r).day(d).begin(b).powerBegin');
                    power_freq = group(g).rat(r).day(d).begin(b).power_freq;
                    [ci, ~] = bootci(2000, {@mean, samples});

                    % After getting the bootstrap samples, plot them
                    lower = ci(1,:);
                    upper = ci(2,:);
                                        
                    patch([power_freq' fliplr(power_freq')], [lower fliplr(upper)], colors(g), 'EdgeColor','none', 'FaceAlpha', 0.3);
                    plot(power_freq, (upper+lower)/2, colors(g))
                    axis square; xlim([0, 45]); xline(2, '--');
                    xlabel('Frequency [Hz]')
                    ylabel('Power [dB]')
                end
                legend({'apoE4','', '','','','','','','','','','', 'apoE3'})
            end
        end
    end
    %saveas(gcf, [saveDir filesep '20230602_PSD.png'])
end