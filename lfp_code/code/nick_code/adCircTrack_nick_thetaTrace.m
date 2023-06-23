function group = adCircTrack_nick_thetaTrace(group)
    saveDir = 'C:\Users\nrive\Projects\Colgin Lab\apoE4\figures\lfp_e3e4';
    colors = 'bk';
    for g = 1:length(group)
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                figure(d);
                if(g == 1)
                    clf;
                end

                for b = 1:length(group(g).rat(r).day(d).begin)
                    thetaTime = group(g).rat(r).day(d).begin(b).thetaTime;
                    for c = 1:length(group(g).rat(r).day(d).begin(b).lapCuts)
                        thetaTrace = group(g).rat(r).day(d).begin(b).lapCuts(c).thetaCuts;

                        subplot(4,1,b); hold on; 
                        
                     % Also get spec per lap
%                         Pxx = get_wavelet_power(mean(thetaTrace,2), 2000, [1, 25], 6);
%                         if(g == 1)
%                             offset = -0.35;
%                         else
%                             offset = 0.35;
%                         end
%                         colormap(hot);
%                         imagesc(thetaTime+((2*c-1)+offset), 1:50, pow2db(Pxx)); colorbar; title(['apoE4(L) and apoE3(R) per lap Day ' num2str(d) ' Begin ' num2str(b)])
%                         ylabel('Frequency')
%                         xlabel('Lap # [-0.3, +0.3] window from onset')
%                         xticks(1:2:length(group(g).rat(r).day(d).begin(b).lapCuts)*2)
%                         xticklabels(1:length(group(g).rat(r).day(d).begin(b).lapCuts))
                        
                        plot(thetaTime+(c), mean(thetaTrace,2), colors(g))
                        title(['apoE4 (Blue) and apoE3 (Black) Waveform per lap Day ' num2str(d) ' Begin ' num2str(b)])
                        xlabel('Lap# + Time from theta trough [s]')
                        ylabel('Z-scored Amplitude')
                    end
                end
                %saveas(gcf, [saveDir filesep '20230615_waveSpec_Day' num2str(d) '.png'])
                saveas(gcf, [saveDir filesep '20230615_rawTracePerLap_Day' num2str(d) '.png'])
            end
        end
    end
end