for g = 1:length(rat381LightsOffrat326z)
    for r = 1:length(rat381LightsOffrat326z(g).rat)
        for d = 1:length(rat381LightsOffrat326z(g).rat(r).day)
            for b = 1:length(rat381LightsOffrat326z(g).rat(r).day(d).begin)
                coords = rat381LightsOffrat326z(g).rat(r).day(d).begin(b).coords;

                subplot(2,1,1);
                plot(coords(:,2), coords(:,3))
                axis square;

                coords = remove_artifacts_nick(coords(:,2:3), 5);
                subplot(2,1,2);
                plot(coords(:,1), coords(:,2))
                axis square;
            end
        end
    end
end