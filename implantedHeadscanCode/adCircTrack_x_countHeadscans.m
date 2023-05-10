figure(1); clf;
figure(2); clf;
figure(3); clf;
FPS = 30;

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        headscan_count = zeros(length(group(g).rat(r).day), 4);
        headscan_duration = zeros(length(group(g).rat(r).day), 4);
        headscan_time = zeros(length(group(g).rat(r).day), 4);

        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            for b = 1:length(group(g).rat(r).day(d).begin)
                pss = group(g).rat(r).day(d).begin(b).pss;
                filt_pss = group(g).rat(r).day(d).begin(b).filt_pss;

                headscan_count(d,b) = group(g).rat(r).day(d).begin(b).num_headscans;
                headscan_duration(d,b) = mean(filt_pss(:,2)-filt_pss(:,1));
                headscan_time(d,b) = sum(pss)/FPS;
            end
        end

        figure(1); hold on;
        x = (1:length(group(g).rat(r).day))+(0.1*g - 0.05);
        y = mean(headscan_count,2);
        err = std(headscan_count, [], 2) / sqrt(length(group(g).rat(r).day));
        errorbar(x, y, err)

        figure(2);  hold on;
        x = (1:length(group(g).rat(r).day))+(0.1*g - 0.05);
        y = mean(headscan_duration,2);
        err = std(headscan_duration, [], 2) / sqrt(length(group(g).rat(r).day));
        errorbar(x, y, err)

        figure(3); hold on;
        x = (1:length(group(g).rat(r).day))+(0.1*g - 0.05);
        y = mean(headscan_time,2);
        err = std(headscan_time, [], 2) / sqrt(length(group(g).rat(r).day));
        errorbar(x, y, err)
    end
end

figure(1);
title('Number of Headscans')
xlabel('Day');
ylabel('Number of headscans [# of events]');
xticks(1:5)
legend({'apoE3', 'apoE4'})

figure(2);
title('Average Headscan Duration')
xlabel('Day');
ylabel('Duration [s]');
xticks(1:5)
legend({'apoE3', 'apoE4'})

figure(3);
title('Total Headscan Time')
xlabel('Day');
ylabel('Time spent headscanning per begin [s]');
xticks(1:5)
legend({'apoE3', 'apoE4'})