tic
opa_build_struct;

%% During session 1, compare exploration times between novel and familiar
% (Even though this is a control)
session_number = 'S2';
event_to_track = 'first3_objA';
session_spec = opa_struct(strcmp({opa_struct.sessionNumber}, session_number));

% Get the values per animal, per minute, per type, per group
max_animals = 6;
max_minutes = 1;
session_types = unique({session_spec.sessionType});
groups = unique({session_spec.group});
s1_struct = {};
s1_index = 1;

for g = 1:length(groups)
    for t = 1:length(session_types)
        temp_session_spec = session_spec(strcmp({session_spec.group}, groups{g}));
        s1_struct(s1_index).group_name = groups{g};
        s1_struct(s1_index).type = session_types{t};

        temp_session_spec = temp_session_spec(strcmp({temp_session_spec.sessionType}, s1_struct(s1_index).type));
        % Get all values for a single group/session type for plotting
        objA_struct = nan(max_animals, max_minutes);
        for rat = 1:length(temp_session_spec)
            %objA_struct(rat, 1:temp_session_spec(rat).vidMinutes) = temp_session_spec(rat).(event_to_track);
            objA_struct(rat) = temp_session_spec(rat).(event_to_track);
        end
        s1_struct(s1_index).mean = mean(objA_struct, 'omitnan');
        s1_struct(s1_index).sem = std(objA_struct, 'omitnan') / sqrt(length(temp_session_spec));
        s1_index = s1_index+1;
    end
end

%% Plot the different object pair types
for i = 1:6
    figure(i); clf; hold on;
    e3 = errorbar(s1_struct(i).mean, s1_struct(i).sem);
    e4 = errorbar(s1_struct(i+6).mean, s1_struct(i+6).sem);
    e3.Color = 'k';
    e4.Color = 'b';
    title_string = [session_types{i} ' Session ' session_number];
    file_string = strrep(title_string, ' ', '_');

    title(title_string)
    ylabel('Number of Interactions')
    xlabel('Time in session [min]')
    legend({'apoE3', 'apoE4'})
    %saveas(gcf, ['../figures/numEventsObjB/20230407_' file_string '.png'])
end
toc