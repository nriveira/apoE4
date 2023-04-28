opa_build_struct;

%% During session, compare exploration times to practice session with no objects
% First get practice session
 session_number = 'S3';
practice_number = 'P3';
firstXframes = 1800 * 3; % 30 fps -> 1800 frames 60 seconds

% Isolate sessions and practice sessions
session_prac = opa_struct(strcmp({opa_struct.sessionNumber}, practice_number));
session_spec = opa_struct(strcmp({opa_struct.sessionNumber}, session_number));

% Get the values per animal, per minute, per type, per group
max_animals = 6;
session_types = unique({session_spec.sessionType});
groups = unique({session_spec.group});
s_struct = {};
s_index = 1;

for g = 1:length(groups)
    for t = 1:length(session_types)
        group_session_spec = session_spec(strcmp({session_spec.group}, groups{g}));
        s_struct(s_index).group_name = groups{g};
        s_struct(s_index).type = session_types{t};

        temp_session_spec = group_session_spec(strcmp({group_session_spec.sessionType}, s_struct(s_index).type));

        % Get all values for a single group/session type for plotting
        obj_struct = nan(length(temp_session_spec), 3);
        for rat = 1:length(temp_session_spec)
            temp_session_prac = session_prac(strcmp({session_prac.name}, temp_session_spec(rat).name));
            if(~isempty(temp_session_prac))
                objA_coords = [temp_session_spec(rat).session.objA_x_coord temp_session_spec(rat).session.objA_y_coord];
                objB_coords = [temp_session_spec(rat).session.objB_x_coord temp_session_spec(rat).session.objB_y_coord];
    
                distance_objA = sqrt((temp_session_prac.session.nose_x_coord-objA_coords(1)).^2 + (temp_session_prac.session.nose_y_coord-objA_coords(2)).^2)./parameters.est_conversion;
                distance_objB = sqrt((temp_session_prac.session.nose_x_coord-objB_coords(1)).^2 + (temp_session_prac.session.nose_y_coord-objB_coords(2)).^2)./parameters.est_conversion;
    
                prac_time_objA = sum(distance_objA(1:firstXframes) < parameters.cm_from_object)/parameters.fps;
                prac_time_objB = sum(distance_objB(1:firstXframes) < parameters.cm_from_object)/parameters.fps;
    
                sess_time_objA = sum(temp_session_spec(rat).perMin_objA(1:3));
                sess_time_objB = sum(temp_session_spec(rat).perMin_objB(1:3));                
    
                obj_struct(rat,1) = sess_time_objA / (prac_time_objA+sess_time_objA); % Gives time in seconds
                obj_struct(rat,2) = sess_time_objB / (prac_time_objB+sess_time_objB);
                obj_struct(rat,3) = sess_time_objA / (sess_time_objA+sess_time_objB);
            end
        end
        s_struct(s_index).obj_struct = obj_struct;
        s_struct(s_index).mean = mean(obj_struct, 'omitnan');
        s_struct(s_index).sem = std(obj_struct, 'omitnan') ./ sqrt(length(temp_session_spec));
        s_index = s_index+1;
    end
end

%% Plot Data
figure(1); clf;
data = reshape([s_struct(:).mean], [3,12]); 
bar(data);
xticklabels({'Object A','Object B','Obj A-B'})
ylabel('Discrimination index')
legend(session_types)

