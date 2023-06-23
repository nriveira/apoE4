function dlc_struct = get_dlc_struct(inputFile)
% get_dlc_struct Build a matlab compatable stuct for the output of DLC
%   Output is DLC file with rectangular coordinates and the predicted
%   center (for circle track)
    dlc_output = readmatrix(inputFile);
    dlc_no_scorer = dlc_output(:,2:end);
    pctile = 0.99;

    dlc_struct = [];
    % Add each body part to a different row of DLC struct, points are
    % centered
    for i = 1:size(dlc_no_scorer, 2)/3
        dlc_struct(i).bodypart = 'Unnamed';
        x_coord = dlc_no_scorer(:,(i-1)*3+1);
        y_coord = dlc_no_scorer(:,(i-1)*3+2);
        dlc_struct(i).probability = dlc_no_scorer(:,(i-1)*3+3);

        % Use high confidence points to center track
        conf = dlc_struct(i).probability > prctile(dlc_struct(i).probability, pctile);
        x_conf = x_coord(conf, :);
        y_conf = y_coord(conf, :);

        % Convert to polar coordinates after centering
        dlc_struct(i).predicted_center_x = ((max(x_conf)+min(x_conf))/2);
        dlc_struct(i).predicted_center_y = ((max(y_conf)+min(y_conf))/2);

        dlc_struct(i).x_centered = x_coord - dlc_struct(i).predicted_center_x;
        dlc_struct(i).y_centered = y_coord - dlc_struct(i).predicted_center_y;

        distance_traveled = (sqrt(diff(dlc_struct(i).x_centered).^2 + diff(dlc_struct(i).y_centered).^2))/30;
        prev_x = dlc_struct(i).x_centered(1);
        prev_y = dlc_struct(i).y_centered(1);
 
        % Smooth the data
        for a = 2:length(dlc_struct(i).x_centered)
            if(distance_traveled(a-1) > 20) %~25 cm/s jump in one frame
                dlc_struct(i).x_centered(a) = prev_x;
                dlc_struct(i).y_centered(a) = prev_y;
            end
            prev_x = dlc_struct(i).x_centered(a);
            prev_y = dlc_struct(i).y_centered(a);
        end    

        dlc_struct(i).dadt = (sqrt(diff(dlc_struct(i).x_centered).^2 + diff(dlc_struct(i).y_centered).^2))/30;
        dlc_struct(i).polar_coord = pi+atan2(dlc_struct(i).y_centered, dlc_struct(i).x_centered);
        dlc_struct(i).predicted_r = sqrt(dlc_struct(i).x_centered.^2 + dlc_struct(i).y_centered.^2);

        % Centering first, by finding the median radius along 32 split bins
        [~, ~, bin] = histcounts(dlc_struct(i).polar_coord, 31);
        bin_centers = zeros(32,1);
        for bins = 1:length(bin_centers)
            bin_centers(bins) = median(dlc_struct(i).predicted_r(bin(bin==bins)));
        end
        dlc_struct(i).r_length = mean(bin_centers);
    end
end