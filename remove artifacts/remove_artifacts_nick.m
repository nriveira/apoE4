function xy_coordinates = remove_artifacts_nick(xy_coordinates, threshold)
%REMOVE_ARTIFACTS_NICK Based on the velocity traveled between frames,
%remove jump discontinuities
    distance_traveled = sqrt(diff(xy_coordinates(:,1)).^2 + diff(xy_coordinates(:,2)).^2);
    prev = xy_coordinates(1,:);

    for i = 2:length(xy_coordinates)
        if(distance_traveled(i-1) > threshold)
            xy_coordinates(i-1,:) = prev;
        end
        prev = xy_coordinates(i-1,:);
    end
end

