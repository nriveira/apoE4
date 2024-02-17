fileLoc = dir(fullfile('../headscanData/centerfit_videos/', '**/*.mp4'));
center.vals = zeros(length(fileLoc)-2, 2);
center.filenames = {};

%% Create center fit values and names
for f = 1:length(fileLoc)
    videoDetails = extractAfter(fileLoc(f).folder, 'centerfit_videos\');
    [center, radius] = fitcenter([fileLoc(f).folder filesep fileLoc(f).name]);
    center.vals(f-2, :) = center;
    center.filenames{f-2} = videoDetails;
end

%% Format strings to match lookup from other program
center_lookup = {};

for f = 1:length(center.file_names)
    spltStr = split(center.file_names{f}, '\');
    if(length(spltStr) == 4)
        formString = [spltStr{2} '_d' extractAfter(spltStr{3}, 'day') spltStr{4}];
        center_lookup(f).str = formString;
        center_lookup(f).center = center.center_vals(f,:);
    else
        formString = [extractBefore(spltStr{2}, '_') '_' spltStr{3}];
        center_lookup(f).str = formString;
        center_lookup(f).center = center.center_vals(f,:);
    end
end

%% Manually fix the last ones (The unnammed ones)
center_lookup(107).str = 'rat381_d1s1';
center_lookup(108).str = 'rat381_d1s2';
center_lookup(109).str = 'rat381_d1s3';
center_lookup(110).str = 'rat381_d2s1';
center_lookup(111).str = 'rat381_d2s2';
center_lookup(112).str = 'rat381_d2s3';
center_lookup(113).str = 'rat381_d3s1';
center_lookup(114).str = 'rat381_d3s2';
center_lookup(115).str = 'rat381_d3s3';
center_lookup(116).str = 'rat381_d4s1';
center_lookup(117).str = 'rat381_d4s2';
center_lookup(118).str = 'rat381_d4s3';
center_lookup(119).str = 'rat382_d1s1';
center_lookup(120).str = 'rat382_d1s2';
center_lookup(121).str = 'rat382_d1s3';
center_lookup(122).str = 'rat382_d2s1';
center_lookup(123).str = 'rat382_d2s2';
center_lookup(124).str = 'rat382_d2s3';
center_lookup(125).str = 'rat382_d3s1';
center_lookup(126).str = 'rat382_d3s2';
center_lookup(127).str = 'rat382_d3s3';
center_lookup(128).str = 'rat382_d4s1';
center_lookup(129).str = 'rat382_d4s2';
center_lookup(130).str = 'rat382_d4s3';