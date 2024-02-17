function dlc_struct = get_dlc_struct(inputFile)
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

        load('center_lookup.mat');
        % Use the lookup index to get the center values for that video
        lookup_val = extractBetween(inputFile, 'iteration1\', 'DLC_resnet50');
        if(sum(strcmp({center_lookup.str}, lookup_val)) == 1)
            center_val = center_lookup(strcmp({center_lookup.str}, lookup_val)).center;
            xc = center_val(1);
            yc = center_val(2);
        else 
            lookup_val
        end

        % Find the center to center all xy values
        [xc_old,yc_old,R,~] = circfit(x_coord, y_coord);
        [abs(xc - xc_old) abs(yc-yc_old)];

        % Convert to polar coordinates after centering
        dlc_struct(i).predicted_center_x = xc;
        dlc_struct(i).predicted_center_y = yc;

        x_centered = x_coord - xc;
        y_centered = y_coord - yc;
        xy_centered = [x_centered y_centered];
        xy_centered = remove_artifacts_nick(xy_centered);

        dlc_struct(i).x_centered = xy_centered(:,1);
        dlc_struct(i).y_centered = xy_centered(:,2);
       
        [polar_coord, dlc_struct(i).predicted_r] = cart2pol(dlc_struct(i).x_centered, dlc_struct(i).y_centered);
        dlc_struct(i).polar_coord = polar_coord;
        dlc_struct(i).dadt = diff(unwrap(dlc_struct(i).polar_coord)); % Gives dadt in pixels / frame
        dlc_struct(i).r_length = R;
    end
end