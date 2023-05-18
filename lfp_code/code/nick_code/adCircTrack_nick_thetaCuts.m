function group = adCircTrack_nick_thetaCuts(group)
    figure(1); clf;
    for g = 1:length(group)
        subplot(2,1,g); hold on;
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                samples = [];
                for b = 1:length(group(g).rat(r).day(d).begin)  
                    samples = [samples, group(g).rat(r).day(d).begin(b).thetaCuts];
                end
                plot(group(g).rat(r).day(d).begin(b).thetaTime, mean(samples,2))
            end
        end
    end
end