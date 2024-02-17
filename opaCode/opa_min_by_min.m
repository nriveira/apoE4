%% Analyze session struct
apoE3 = {'rat385','rat386','rat387','rat388', 'rat389', 'rat390'};
analyzed_struct = opa_analyze_session_struct(opa_struct, apoE3);

% Cycle through all four groups and
% Calculate the mean and standard error of all data
figure_data = {};
figure_data.apoE3 = get_fig_values(analyzed_struct.apoE3_struct, opa_struct);
figure_data.apoE4 = get_fig_values(analyzed_struct.apoE4_struct, opa_struct);

groups = {'Familiar', 'NO', 'NL', 'NO NL'};
figure(1); clf; hold on;
for sub = 1:4
    subplot(4,1,sub); 
    bar([figure_data.apoE3(sub).perMin_mean, figure_data.apoE4(sub).perMin_mean])
    legend({'apoE3','apoE4'}, 'Location','bestoutside')
    title(groups{sub});
    ylabel('Time near object [s]')
    xlabel('Minutes of session')
    xlim([0, 10])
end

figure(2); clf; hold on;
for sub = 1:4
    subplot(4,1,sub); 
    bar([figure_data.apoE3(sub).numEvents_mean, figure_data.apoE4(sub).numEvents_mean])
    legend({'apoE3','apoE4'}, 'Location','bestoutside')
    title(groups{sub});
    ylabel('Number of events')
    xlabel('Minutes of session')
    xlim([0, 10])
end

function apo = get_fig_values(analyzed_struct_apo, opa_struct)
    apo = [];
    for i = 1:4
        apo(i).name = analyzed_struct_apo(i).group;
        perMin_data = nan(max([opa_struct.vidMinutes]), length(analyzed_struct_apo(i).struct));
        event_data = nan(max([opa_struct.vidMinutes]), length(analyzed_struct_apo(i).struct));
    
        % Cycle through each session and pull out the events
        % Tricky because they are not all the same length videos!
        for s = 1:length(analyzed_struct_apo(i).struct)
            perMin_s = analyzed_struct_apo(i).struct(s).perMin_objA;
            numEvents_s = analyzed_struct_apo(i).struct(s).numEvents_objA;
    
            perMin_data(1:length(perMin_s), s) = perMin_s;
            event_data(1:length(perMin_s), s) = numEvents_s;
        end

        apo(i).perMin_data = perMin_data;
        apo(i).numEvents_data = event_data;
    
        apo(i).perMin_mean = mean(perMin_data, 2);
        apo(i).numEvents_mean = mean(event_data, 2);
    
        apo(i).perMin_sem = sqrt(var(perMin_data, 0, 2)./size(perMin_data,2));
        apo(i).numEvents_sem = sqrt(var(event_data, 0, 2)./size(event_data,2));
    end
end