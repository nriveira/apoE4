fileLoc = '../headscanData/iteration1';
fileDir = dir(fileLoc);

headscan_struct = {};
apoE4 = {'rat377', 'rat378', 'rat379', 'rat381', 'rat382'};

%% Load data once, then use the structure
for file = 3:length(fileDir)
    fileName = [fileDir(file).folder filesep fileDir(file).name];
    headscan = analyze_headscans(fileName, param);

    headscan_struct(file-2).name = extractBefore(fileDir(file).name, '_d');
    headscan_struct(file-2).day = extractBetween(fileDir(file).name, '_d', 's');
    headscan_struct(file-2).session = extractBetween(fileDir(file).name, 's', 'DLC_resnet50');
    if any(strcmp(apoE4, headscan_struct(file-2).name))
        headscan_struct(file-2).group = 'apoE4';
    else
        headscan_struct(file-2).group = 'apoE3';
    end
    
    headscan_struct(file-2).headscan = headscan;
end

%% Using data structure, split all of the data being plotted to get lap info
apoE3_headscans = headscan_struct(strcmp({headscan_struct(:).group}, 'apoE3'));
apoE4_headscans = headscan_struct(strcmp({headscan_struct(:).group}, 'apoE4'));

% Aggregated values to plot
e3_vals = zeros(12, 6, 10);
e4_vals = zeros(12, 5, 10);

% Values represent
% - Laps per session
% - Time per lap (Total time / Number of laps)
% - Pause time per lap (Total pause time / Number of laps)
% - Headscans per lap (Number of headscan events)
% - Pause time per lap (normalized by pause time)
% - Headscan time per lap
% - Average Headscan time
% - Headscan ISI (Total time / number of scans)
% - Scan rate (Head scan events / pause time)

for day = 1:4
    % Split the days back up
    dayE3 = apoE3_headscans(strcmp([apoE3_headscans(:).day], num2str(day)));
    dayE4 = apoE4_headscans(strcmp([apoE4_headscans(:).day], num2str(day)));

    for sess = 1:3
        % Also gather by sessions
        sessionE3 = dayE3(strcmp([dayE3(:).session], num2str(sess)));
        sessionE4 = dayE4(strcmp([dayE4(:).session], num2str(sess)));
        index = 3*(day-1)+sess;
        
        % For each group, gather measures to plot
        for apoE3 = 1:length(sessionE3)
            variablesE3 = get_figure_data(sessionE3(apoE3));
            e3_vals(index, apoE3, :) = variablesE3;
        end
        for apoE4 = 1:length(sessionE4)
            variablesE4 = get_figure_data(sessionE4(apoE4));
            e4_vals(index, apoE4, :) = variablesE4;
        end
    end
end

%% Plot all values
titles = {'Laps Ran','Total time per Number of laps','Pause time per Number of laps','Number of Headscans per Number of laps','Pause time per Total time','Percent Time Headscanning', 'Average Headscan Duration', 'Headscan Inter-scan Interval Time', 'Headscan Rate', 'Number of Headscans'};
ylab = {'Laps','Time [s per lap]','Pause time [s per lap]', 'Number of headscans per lap', 'Pause time [%]', 'Headscan time [%]', 'Average headscan duration [s]', 'Inter-scan Interval [s]', 'Headscans/Pause Time', 'Headscans (Count)'};
figureSave = '../figures/poster figures';

for feature = 1:size(e3_vals, 3)
    figure(feature); clf; hold on;
    
    e3mean = mean(e3_vals(:,:,feature), 2);
    e3error = std(e3_vals(:,:,feature), 0, 2)/sqrt(size(e3_vals(:,:,feature), 2));
    
    e4mean = mean(e4_vals(:,:,feature), 2);
    e4error = std(e4_vals(:,:,feature), 0, 2)/sqrt(size(e4_vals(:,:,feature), 2));
    
    errorbar(e3mean, e3error, 'k');
    errorbar(e4mean, e4error, 'b');
    
    plot(e3_vals(:,:,feature),'.k');
    plot(e4_vals(:,:,feature), '.b');
    
    title(titles{feature})
    xticks(0:13)
    xlim([0,13])
    xticklabels({'','Day 1','S2','S3','Day 2','S2','S3','Day 3','S2','S3','Day 4','S2','S3'})
    xlabel('Sessions')
    ylabel(ylab{feature})
    ylim([0 inf])
    legend({'apoE3', 'apoE4'})
    saveas(gcf, [figureSave filesep '20230623_' char(strrep(titles(feature), ' ', '_'))], 'png')
end

% Values represent
% - Laps per session
% - Time per lap (Total time / Number of laps)
% - Pause time per lap (Total pause time / Number of laps)
% - Headscans per lap (Number of headscan events)
% - Pause time percentage (Pause time / Total Time)
% - Headscan time (Total time spent headscanning)
% - Average Headscan time (Average headscan duration)
% - Headscan ISI (Total time / number of scans)
% - Scan rate (Head scan events / pause time)
% - Number of events (Total number of headscan events)

function variables = get_figure_data(session)
    num_laps = session.headscan.num_laps+1;
    time_per_lap = session.headscan.total_time / (num_laps);
    pause_per_lap = session.headscan.pause_time / (num_laps);
    headscans_per_lap = session.headscan.num_headscans / (num_laps);

    pause_time_percent = session.headscan.pause_time / session.headscan.total_time;
    % Calculate the amount of time spent headscanning (remember fps)
    headscan_time = sum(session.headscan.filt_pss(:,2)-session.headscan.filt_pss(:,1))/session.headscan.total_time;
    avg_headscan_time = mean(session.headscan.filt_pss(:,2)-session.headscan.filt_pss(:,1))/30;
    headscan_isi = session.headscan.total_time / session.headscan.num_headscans;
    headscan_rate = session.headscan.num_headscans / session.headscan.pause_time;
    num_events = session.headscan.num_headscans;

    variables = [num_laps-1, ...
                 time_per_lap, ...
                 pause_per_lap, ...
                 headscans_per_lap, ...
                 pause_time_percent, ... 
                 headscan_time, ...
                 avg_headscan_time, ...
                 headscan_isi, ...
                 headscan_rate, ...
                 num_events];
end