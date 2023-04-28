function headscans = dlc_get_runspeeds(dataLoc, param)
    headscan_dir = dir(dataLoc);
    headscans = {};
    
    for i = 3:length(headscan_dir)
        fileLoc = [headscan_dir(i).folder filesep headscan_dir(i).name];
        headscan_struct = analyze_headscans(fileLoc, param);
    
        headscans(i-2).dlc_struct = headscan_struct.dlc_struct;
        headscans(i-2).dadt = headscan_struct.dadt;
        headscans(i-2).polr = headscan_struct.polar_r_in_cm;
        headscans(i-2).drdt = headscan_struct.drdt;
        headscans(i-2).pss = headscan_struct.pss;
    
        headscans(i-2).laps = headscan_struct.laps;
        headscans(i-2).headscan_time = sum(headscans(i-2).filt_pss)/param.fps;
    end
end