function xy_coordinates = remove_artifacts_nick(xy_coordinates)
%REMOVE_ARTIFACTS_NICK Based on the velocity traveled between frames,
%remove jump discontinuities
                       % sqrt((x2-x1)^2 + (y2-y1)^2)
    distance_traveled = sqrt(diff(xy_coordinates(:,1)).^2 + diff(xy_coordinates(:,2)).^2);
    threshold = prctile(distance_traveled, 95);

    xy_coordinates(distance_traveled > threshold, :) = nan;
    %xy_coordinates = movmedian(xy_coordinates, 15, 1, 'omitnan');

    prev = xy_coordinates(1,:);
    for i = 2:length(xy_coordinates)
        if(isnan(xy_coordinates(i,:)))
            xy_coordinates(i,:) = prev;
        end
        prev = xy_coordinates(i,:);
    end
end