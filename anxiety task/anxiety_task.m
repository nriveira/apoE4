%load('anxiety data.mat');

groups = anxiety_test(:,1);
roi = anxiety_test(:,2);
time_spent = anxiety_test(:,3);
times_crossed = anxiety_test(:,4);

figure(1); clf; hold on;
bar_means = ones(5,2);
bar_sterr = ones(5,2);

for g = 1:2
    for r = 1:5
        samples = times_crossed([groups == g] & [roi == r]);
        bar_means(r,g) = mean(samples);
        bar_sterr(r,g) = std(samples)/sqrt(length(samples));
    end
end

x = 1:5;
bar(bar_means);
errorbar(x-0.15, bar_means(:,1), bar_sterr(:,1), 'k', 'LineStyle','none');
errorbar(x+0.15, bar_means(:,2), bar_sterr(:,2), 'k', 'LineStyle','none');

for g = 1:2
    for r = 1:5
        samples = times_crossed([groups == g] & [roi == r]);
        if(g == 1)
            scatter(r*ones(size(samples))-0.15, samples, '.k')
        else
            scatter(r*ones(size(samples))+0.15, samples, '.k')
        end
    end
end

title('Times Crossed ROI Boundary');
ylabel('Number of Crosses');
xlabel('Region of Interest');
xticklabels({'Center', '', 'Top', '', 'Bottom', '', 'Left', '', 'Right'})

legend('apoE3', 'apoE4')