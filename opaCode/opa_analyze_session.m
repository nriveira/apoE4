function session = opa_analyze_session(file, parameters, plot_session)
    dlc_struct = opa_get_dlc_struct(file, parameters.bodyparts);
    
    % Plot trajectory of rat during this session
    session.nose_x_coord = medfilt1(dlc_struct(strcmp({dlc_struct(:).bodypart},'nose')).x_coord);
    session.nose_y_coord = medfilt1(dlc_struct(strcmp({dlc_struct(:).bodypart},'nose')).y_coord);
    
    % Estimate object A and B positions    
    session.objA_x_coord = median(dlc_struct(strcmp({dlc_struct(:).bodypart},'objectA')).x_coord);
    session.objA_y_coord = median(dlc_struct(strcmp({dlc_struct(:).bodypart},'objectA')).y_coord);
    
    session.objB_x_coord = median(dlc_struct(strcmp({dlc_struct(:).bodypart},'objectB')).x_coord);
    session.objB_y_coord = median(dlc_struct(strcmp({dlc_struct(:).bodypart},'objectB')).y_coord);
    
    % Find difference between points
    session.nose_objA_dist = sqrt((session.nose_x_coord-session.objA_x_coord).^2 + (session.nose_y_coord-session.objA_y_coord).^2)./parameters.est_conversion;
    session.nose_objB_dist = sqrt((session.nose_x_coord-session.objB_x_coord).^2 + (session.nose_y_coord-session.objB_y_coord).^2)./parameters.est_conversion;
    
    if(plot_session)
        plot_functions(session, parameters);
    end
end

% Helper function that plots all of the results for a single
% session
function output = plot_functions(session, parameters)  
    % Create vis of rat trajectory and object locations
    plot(session.nose_x_coord, session.nose_y_coord, '.k');

    plot(session.objA_x_coord, session.objA_y_coord, 'xb')
    plot(session.objB_x_coord, session.objB_y_coord, 'xg')
    
    % Draw area considered 'object exploration events'
    viscircles([session.objA_x_coord, session.objA_y_coord], parameters.cm_from_object*parameters.est_conversion);
    viscircles([session.objB_x_coord, session.objB_y_coord], parameters.cm_from_object*parameters.est_conversion);

    plot(session.nose_x_coord(session.nose_objA_dist < parameters.cm_from_object), session.nose_y_coord(session.nose_objA_dist < parameters.cm_from_object), '.b')
    plot(session.nose_x_coord(session.nose_objB_dist < parameters.cm_from_object), session.nose_y_coord(session.nose_objB_dist < parameters.cm_from_object), '.g')
end