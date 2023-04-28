clear all; close all; clc;
% load('I:\JBT_RERUN\DATASET\Rerun_Uber_Data_Struct_with3cmBinRm_withStateTimes_withSTCorrs_withSpatMets_09062017.mat')
% cellRegion = corrProj_6_reduceMetricsForEachCellPair(region);
load('I:\JBT_RERUN\DATASET\Rerun_CellPair_Data_Struct_with3cmBins_092017.mat')
lags = -5000:2:5000;
for ii=1:226
    run_xcorr = cellRegion(1).cellPair(ii).state(1).normStXCorr{1};
    rem_xcorr = cellRegion(1).cellPair(ii).state(2).normStXCorr{1};
    cell_pair_indexer = cellRegion(1).cellPair(ii).info;
    rsp_mag = cellRegion(1).cellPair(ii).state(1).relSpatPhiMag / sqrt(2);
    
    % Animal ID & Session ID
    switch cell_pair_indexer(1)
        case 1
            animal_id = 'Rat_109';
            switch cell_pair_indexer(2)
                case 1
                    session_id = '7/30';
                case 2
                    session_id = '8/2';
                case 3
                    session_id = '8/5';
                case 4
                    session_id = '8/17';
                case 5
                    session_id = '8/21';
                case 6
                    session_id = '8/30';
                case 7
                    session_id = '9/2';
            end
        case 2
            animal_id = 'Rat_20';
            switch cell_pair_indexer(2)
                case 1
                    session_id = '12/7';
                case 2
                    session_id = '12/9';
                case 3
                    session_id = '12/13';
                case 4
                    session_id = '12/15';
                case 5
                    session_id = '12/17';
            end
        case 3
            animal_id = 'Rat_26';
            switch cell_pair_indexer(2)
                case 1
                    session_id = '8/5';
                case 2
                    session_id = '8/8';
                case 3
                    session_id = '8/12';
                case 4
                    session_id = '8/14';
            end
        case 4
            animal_id = 'Rat_29';
            session_id = '10/24';
        case 5
            animal_id = 'Rat_74';
            switch cell_pair_indexer(2)
                case 1
                    session_id = '4/2';
                case 2
                    session_id = '4/6';
                case 3
                    session_id = '4/11';
                case 4
                    session_id = '4/14';
                case 5
                    session_id = '4/18';
                case 6
                    session_id = '4/26';
                case 7
                    session_id = '4/30';
            end
        case 6
            animal_id = 'Rat_78';
            switch cell_pair_indexer(2)
                case 1
                    session_id = '7/24';
                case 2
                    session_id = '8/5';
                case 3
                    session_id = '8/9';
                case 4
                    session_id = '8/12';
                case 5
                    session_id = '8/17';
            end
    end
    
    % Cell IDs
    cell_1_id = ['TT_',num2str(cell_pair_indexer(4)),'_',num2str(cell_pair_indexer(5)),'.t'];
    cell_2_id = ['TT_',num2str(cell_pair_indexer(6)),'_',num2str(cell_pair_indexer(7)),'.t'];
    
    %Plot and save
    ffa = figure('Units','inches','Position',[0 0 11 8.5], 'PaperSize',[11 8.5], 'PaperPositionMode','auto', 'PaperUnits','inches'); clf;
        subplot(2,2,1)
        plot(lags, run_xcorr);
        xlim([-250 250])
        title(['RUN -- ',animal_id,', Session: ',session_id,', Cell 1: ',cell_1_id,' vs. Cell 2: ',cell_2_id,'\newlineRelative spatial phase magnitude: ',num2str(rsp_mag,'%1.03f')])
        
        subplot(2,2,3)
        plot(lags, rem_xcorr);
        xlim([-250 250])
        title(['REM -- ',animal_id,', Session: ',session_id,', Cell 1: ',cell_1_id,' vs. Cell 2: ',cell_2_id])
        
        subplot(2,2,2)
        plot(lags, run_xcorr);
        title(['RUN -- ',animal_id,', Session: ',session_id,', Cell 1: ',cell_1_id,' vs. Cell 2: ',cell_2_id,'\newline Relative spatial phase magnitude: ',num2str(rsp_mag,'%1.03f')])
        
        subplot(2,2,4)
        plot(lags, rem_xcorr);
        title(['REM -- ',animal_id,', Session: ',session_id,', Cell 1: ',cell_1_id,' vs. Cell 2: ',cell_2_id])
        
        saveas(ffa,['I:\JBT_RERUN\indiv_cell_plots\RSPM_',num2str(rsp_mag,'%1.03f'),'__Cell_Pair_',num2str(ii),'.png'],'png');
        close all;
end