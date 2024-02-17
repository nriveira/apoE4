function dlc_struct = opa_get_dlc_struct(inputFile, bodyparts)
% get_dlc_struct Build a matlab compatable stuct for the output of DLC
%   Output is DLC file with rectangular coordinates and the predicted
%   center (for circle track)
    opts = detectImportOptions(inputFile);
    dlc_output = readmatrix(inputFile);
    dlc_no_scorer = dlc_output(:,2:end);

    dlc_struct = [];
    % Add each body part to a different row of DLC struct, points are
    % centered
    for i = 1:size(dlc_no_scorer, 2)/3
        dlc_struct(i).bodypart = bodyparts{i};
        xy_coords = dlc_no_scorer(:,(i-1)*3+1:(i-1)*3+2);
        xy_coords = remove_artifacts_nick(xy_coords, 50);
        x_coord = xy_coords(:,1);
        y_coord = xy_coords(:,2);

        dlc_struct(i).x_coord = x_coord;
        dlc_struct(i).y_coord = y_coord;
        dlc_struct(i).probability = dlc_no_scorer(:,(i-1)*3+3);
    end
end