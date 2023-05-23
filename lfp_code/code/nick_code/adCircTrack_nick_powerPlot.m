function group = adCircTrack_nick_powerPlot(group)
    % adCircTrack_nick_powerPlot
    % Plots the bootstrap estimated power across each begin for apoE3 and
    % apoE4 rats
    colors = ['r','g','b','k'];
    saveDir = 'C:\Users\nick\Projects\Colgin Lab\apoE4\figures\lfp_power'; 
    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                % Create new figure per day and give title
                figure; clf; hold on; title([group(g).name ' Day ', num2str(d)])
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
                    axis square; xlim([0, 45]); xline(2, '--');
                    xlabel('Frequency [Hz]')
                    ylabel('Power [dB]')
                end
                legend({'B1','', '','B2','','','B3','','','B4','',''})
                saveas(gcf, [saveDir filesep '20230523_' group(g).name 'Day' num2str(d) '.png'])
            end
        end
    end
end