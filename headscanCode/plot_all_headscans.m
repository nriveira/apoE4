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
e3_vals = zeros(12, 6, 8);
e4_vals = zeros(12, 5, 8);

% Values represent
% - Laps per session
% - Time per lap time
% - Pause time per lap
% - Pause time per lap (Percentage)
% - Headscans per lap (Number of headscan events)
% - Headscan Time (Percentage of pause time)

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
titles = {'Laps Ran','Time per lap','Pause time per lap','Percent time pausing','Headscans per lap','Headscan time','Percent time headscanning', 'Average Headscan Time'};
ylab = {'Laps','Time [s per lap]','Pause time [s per lap]', 'Pause time [percentage]', 'Number of headscans per lap', 'Headscan time [s per lap]', 'Time headscanning [percentage of pauses]', 'Headscan time [s]'};

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
    legend({'apoE3', 'apoE4'})
end

% Values represent
% - Laps per session
% - Time per lap time
% - Pause time per lap
% - Pause time per lap (Percentage)
% - Headscans per lap (Number of headscan events)
% - Headscan Time (Percentage of pause time)
% - Average headscan time
function variables = get_figure_data(session)
    num_laps = session.headscan.num_laps;
    time_per_lap = session.headscan.total_time / (num_laps+1);
    pause_per_lap = session.headscan.pause_time / (num_laps+1);
    pause_per_lap_percentage = session.headscan.pause_time / session.headscan.total_time;
    headscans_per_lap = session.headscan.num_headscans / (num_laps+1);

    % Calculate the amount of time spent headscanning (remember fps)
    headscan_time = sum(session.headscan.filt_pss(:,2)-session.headscan.filt_pss(:,1))/30;
    headscan_time_per_lap = headscan_time / (num_laps+1);
    headscan_time_percentage = headscan_time / session.headscan.total_time;
    avg_headscan_time = mean(session.headscan.filt_pss(:,2)-session.headscan.filt_pss(:,1))/30;

    variables = [num_laps, ...
                 time_per_lap, ...
                 pause_per_lap, ...
                 pause_per_lap_percentage, ... 
                 headscans_per_lap, ...
                 headscan_time_per_lap, ...
                 headscan_time_percentage, ...
                 avg_headscan_time];
end