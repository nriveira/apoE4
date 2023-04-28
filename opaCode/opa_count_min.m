function [time_near_object, num_events] = opa_count_min(distance_array, parameters)
    frames_per_min = parameters.fps*60;
    min_explore_time = parameters.fps * parameters.min_explore_thresh;
    num_minutes = ceil(length(distance_array)/frames_per_min);
    time_near_object = zeros(num_minutes, 1);
    num_events = zeros(num_minutes, 1);
    
    % Filter all exploration times to only include events over threshold
    exploration_distances = [distance_array < parameters.cm_from_object];
    exploration_diff = diff(exploration_distances);

    % Calculate the start and stop time of each event
    starts = find(exploration_diff == 1);
    stops = find(exploration_diff == -1);
    if(length(starts) > length(stops))
        stops(length(starts)) = length(exploration_diff);
    end

    if(length(starts) < length(stops))
        stops = stops(2:end);
    end

    % 
    exploration_events = (stops-starts) < min_explore_time;
    for i = 1:length(exploration_events)
        if(exploration_events(i))
            exploration_distances(starts(i):stops(i)) = 0;
        end
    end

    for i = 1:num_minutes
        frame = ((i-1)*frames_per_min)+1:min(i*frames_per_min, length(exploration_diff));
        time_near_object(i) = sum(exploration_distances(frame))./parameters.fps;
        num_events(i) = sum(exploration_diff(frame) == 1);
    end
end