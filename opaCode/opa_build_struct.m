%% Initialize all parameters
parameters.fps = 30;
parameters.median_filter_length = 8;
parameters.est_conversion = 4.8; % estimated conversion is that the box makes up ~480x480 pixels
parameters.cm_from_object = 15;
parameters.min_explore_thresh = 0.5;
parameters.bodyparts = {'nose','head','objectA','objectB'};
parameters.apoE3 = {'rat385','rat386','rat387','rat388', 'rat389', 'rat390'};

%% Run on all sessions and create a data structure
fileLoc = '../opaData/iteration1';
fileDir = dir(fileLoc);
opa_struct = {};
plot_session = false;

for fileNum = 3:length(fileDir)
    file = [fileLoc filesep fileDir(fileNum).name];
    fileName = extractBefore(fileDir(fileNum).name, 'DLC_resnet50');
    session = opa_analyze_session(file, parameters, plot_session);

    opa_struct(fileNum-2).name = extractBefore(fileName, '_');

    % Logic needed for inputing the groups
    if(any(strcmp(parameters.apoE3, opa_struct(fileNum-2).name)))
        group = 'E3';
    else
        group = 'E4';
    end

    % Assigning all other variables
    opa_struct(fileNum-2).group = group;
    opa_struct(fileNum-2).sessionType = extractAfter(extractBefore(fileName, '_S'), '_');
    opa_struct(fileNum-2).sessionNumber = extractAfter(fileName, [opa_struct(fileNum-2).sessionType, '_']);
    if(isempty(opa_struct(fileNum-2).sessionType))
        opa_struct(fileNum-2).sessionType = opa_struct(fileNum-2).sessionNumber;
    end
    opa_struct(fileNum-2).vidMinutes = [];
    opa_struct(fileNum-2).session = session; % Comment out to reduce data struct size

    [opa_struct(fileNum-2).perMin_objA, opa_struct(fileNum-2).numEvents_objA] = opa_count_min(session.nose_objA_dist, parameters);
    [opa_struct(fileNum-2).perMin_objB, opa_struct(fileNum-2).numEvents_objB] = opa_count_min(session.nose_objB_dist, parameters);

    opa_struct(fileNum-2).objectALoc = [session.objA_x_coord session.objA_y_coord];
    opa_struct(fileNum-2).objectBLoc = [session.objB_x_coord session.objB_y_coord];
    opa_struct(fileNum-2).vidMinutes = length(opa_struct(fileNum-2).numEvents_objB);

    if(opa_struct(fileNum-2).vidMinutes > 3)
        objA_time = sum(opa_struct(fileNum-2).perMin_objA(1:3));
        objB_time = sum(opa_struct(fileNum-2).perMin_objB(1:3));
        % Add in first 3 minutes count to make analyses easier
        opa_struct(fileNum-2).first3_objA = objA_time;
        opa_struct(fileNum-2).first3_objB = objB_time;
        opa_struct(fileNum-2).objA_objB_discrim = objA_time / (objA_time+objB_time);
    end
end