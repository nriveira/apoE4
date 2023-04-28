function analyzed_struct = opa_analyze_session_struct(session_struct, apoGroup)
% ANALYZE_SESSION_STRUCT Build structure for figures
%   Creates an analysis struct based on the individual session analyses
    apoIndex = zeros(length(session_struct), 1);
    for i = 1:length(apoGroup)
        currentRat = apoGroup{i};
        apoIndex(strcmp({session_struct.name}, currentRat)) = 1;
    end

    % Split sessions into E3 and E4 first
    apoE3_struct = session_struct(apoIndex == 1);
    apoE4_struct = session_struct(apoIndex == 0);
    groups = {'NO', 'NL', 'NO_NL'};

    analyzed_struct.apoE3_struct = opa_split_struct(apoE3_struct, groups);
    analyzed_struct.apoE4_struct = opa_split_struct(apoE4_struct, groups);
end

function split_struct = opa_split_struct(apo, groups)
    split_struct = {};
    split_struct(1).group = 'F';
    split_struct(1).struct = apo(contains({apo.sessionType}, 'F'));

    for i = 1:length(groups)
        split_struct(i+1).group = groups{i};
        split_struct(i+1).struct = apo(strcmp({apo.sessionType}, groups{i}));
    end
end