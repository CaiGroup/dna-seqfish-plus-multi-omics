function [offsets] = seqfishtformalign_v3IFext(experimentDir, experimentLabel, ...
    position, numRounds, numCh, ...
    folderArray, chaTform, usechabboffsets, usephysicaloffsets, ...
    refSaveName, hybSaveName, refSaveName2, experimentDir_ref, experimentLabel_ref, hyb1_ref)   
    
    % 2nd phase for decoding points
    %addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
    %addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts\', '-end');
    %addpath('C:\github\streamline-seqFISH\src\beadalignment\', '-end');


    
    %% variables
    chArray = 1:numCh;
    % from DNA seqFISH folder for fiducial and hyb1 offset
    offsetsDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'positions');
    offsetsDir_ref = fullfile(experimentDir_ref, 'analysis', experimentLabel_ref, 'points', 'positions');
    savefnameMat = ['point-offsets_pos' num2str(position) '_'];

    %% get offsets from mat file     
    offsetsName = sprintf(savefnameMat, position);%['20191125_pos' num2str(position) '_offsets_initial_AdjustedPointRatio-logdots-20191213_xyse3_zse3_xyte2_zte2_xyme1_zme1.csv'];
    offsetsPath = getfile(offsetsDir, offsetsName, 'match');
    offsetsPath_ref = getfile(offsetsDir_ref, offsetsName, 'match');

    % get points and apply offsets
    [offsets] = alignpointswrapper_v4IF(chArray, [], offsetsPath, ...
        chaTform, numRounds, folderArray, [], usechabboffsets, usephysicaloffsets, offsetsPath_ref, hyb1_ref); % updated to v3

end