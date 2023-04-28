opa_build_struct;

%% During session, compare exploration times to practice session with no objects
% First get practice session
session_number = 'S1';
practice_number = 'P1';
firstXframes = 1800 * 3; % 30 fps -> 1800 frames 60 seconds

session_prac = opa_struct(strcmp({opa_struct.sessionNumber}, practice_number));
session_spec = opa_struct(strcmp({opa_struct.sessionNumber}, session_number));

% Get the values per animal, per minute, per type, per group
max_animals = 6;
max_minutes = 12;
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
        if(~isempty(temp_session_spec))
            obj_struct = nan(length(temp_session_spec), 1);
            for rat = 1:length(temp_session_spec)
                temp_session_prac = session_prac(strcmp({session_prac.name}, temp_session_spec(rat).name));
                objA_coords = [temp_session_spec(rat).session.objA_x_coord temp_session_spec(rat).session.objA_y_coord];
                objB_coords = [temp_session_spec(rat).session.objB_x_coord temp_session_spec(rat).session.objB_y_coord];
    
                distance = sqrt((temp_session_spec(rat).session.nose_x_coord-objA_coords(1)).^2 + (temp_session_spec(rat).session.nose_y_coord-objA_coords(2)).^2)./parameters.est_conversion;
                prac_time = sum(distance(1:firstXframes) < parameters.cm_from_object)/parameters.fps;
                sess_time = temp_session_spec(rat).perMin_objA(1:3);
    
                obj_struct(rat) = sess_time / (prac_time+sess_time); % Gives time in seconds
            end
            s_struct(s_index).mean = mean(obj_struct, 'omitnan');
            s_struct(s_index).sem = std(obj_struct, 'omitnan') / sqrt(length(temp_session_spec));
            s_index = s_index+1;
        end
    end
end