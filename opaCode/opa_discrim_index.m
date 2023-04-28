opa_build_struct;

%% During session 2, compare exploration times between novel and familiar
% during first 3 minutes
session_number = 'S3';
session_spec = opa_struct(strcmp({opa_struct.sessionNumber}, session_number));

% Get the values per animal, per type, per group
session_types = unique({session_spec.sessionType});
groups = unique({session_spec.group});
s_struct = {};
s_index = 1;

for g = 1:length(groups)
    for t = 1:length(session_types)
        temp_session_spec = session_spec(strcmp({session_spec.group}, groups{g}));
        s_struct(s_index).group_name = groups{g};
        s_struct(s_index).type = session_types{t};

        temp_session_spec = temp_session_spec(strcmp({temp_session_spec.sessionType}, s_struct(s_index).type));
        % Get all values for a single group/session type for plotting
        disc_struct = nan(length(temp_session_spec),1);
        for rat = 1:length(temp_session_spec)
            if(temp_session_spec(rat).vidMinutes > 3)
                time_novel = sum(temp_session_spec(rat).perMin_objB(1:3));
                time_familiar = sum(temp_session_spec(rat).perMin_objA(1:3));
                disc_struct(rat) = time_novel / (time_novel+time_familiar);
            end
        end
        s_struct(s_index).values = disc_struct;
        s_struct(s_index).mean = mean(disc_struct, 'omitnan');
        s_struct(s_index).sem = std(disc_struct, 'omitnan') / sqrt(length(temp_session_spec));
        s_index = s_index+1;
    end
end

figure(1); clf;
bar_array = [s_struct(strcmp({s_struct.group_name}, 'E3')).mean; s_struct(strcmp({s_struct.group_name}, 'E4')).mean];
bar(bar_array')
xticklabels({'F1', 'F2', 'F3', 'NL', 'NO', 'NO+NL'})
ylabel('Discrimination index');
yline(0.5, '--')
legend({'apoE3', 'apoE4'})
title(['Discrimination index during ' session_number])