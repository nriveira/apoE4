function dlc_struct = get_dlc_struct_NHS(inputFile, center, radius)
% get_dlc_struct Build a matlab compatable stuct for the output of DLC
%   Output is DLC file with rectangular coordinates and the predicted
%   center (for circle track)
    dlc_output = readmatrix(inputFile);
    dlc_no_scorer = dlc_output(:,2:end);

    dlc_struct = [];
    % Add each body part to a different row of DLC struct, points are
    % centered
    for i = 1:size(dlc_no_scorer, 2)/3
        dlc_struct(i).bodypart = 'Unnamed';
        x_coord = dlc_no_scorer(:,(i-1)*3+1);
        y_coord = dlc_no_scorer(:,(i-1)*3+2);

        % Use the lookup index to get the center values for that video
        xc = center(1);
        yc = center(2);

        % Find the center to center all xy values
        [xc_old,yc_old,R,~] = circfit(x_coord, y_coord);
        [abs(xc - xc_old) abs(yc-yc_old)];

        x_centered = (x_coord - xc)./(radius/50);
        y_centered = (y_coord - yc)./(radius/50);
        xy_centered = [x_centered y_centered];
        xy_centered = remove_artifacts_nick(xy_centered);
        time = (1:length(x_centered))./30;
        coords = [time', xy_centered];

        dlc_struct(i).x_centered = xy_centered(:,1);
        dlc_struct(i).y_centered = xy_centered(:,2);
        cartesianVelocity = smooth_runspeed(get_runspeed(coords)); % Converting to cm/s

        [polar_coord, predicted_r] = cart2pol(dlc_struct(i).x_centered, dlc_struct(i).y_centered);
        predicted_r(diff(predicted_r) > 100) = nan;
        predicted_r = movmedian(predicted_r, 15, 1, 'omitnan');

        dlc_struct(i).cartVelocity = cartesianVelocity(:,2);
        dlc_struct(i).predicted_r = predicted_r; % Gives predicted r in cm
        dlc_struct(i).polar_coord = pi+polar_coord;
        dlc_struct(i).dadt = diff(unwrap(dlc_struct(i).polar_coord)).*(30*50); % Gives dadt in cm / s
    end
end 