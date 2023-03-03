function offsets = seqfishtformalign_v2lncrna(experimentDir, experimentLabel, ...
    position, numRounds, numCh, ...
    folderArray, chaTform, physicalTform, usechabboffsets, usephysicaloffsets, ...
    refSaveName, hybSaveName, refSaveName2)   
    
    % 2nd phase for decoding points
    %addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
    %addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts\', '-end');
    %addpath('C:\github\streamline-seqFISH\src\beadalignment\', '-end');

    %% variables
    chArray = 1:numCh;
    offsetsDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'positions');
    if exist(offsetsDir, 'dir') ~= 7
        mkdir(offsetsDir);
    end

    savefnameMat = ['point-offsets_pos' num2str(position) '_'];

    % Run the python comand for bead alignment
    
    % new version for more robust alignment correction. 
    % two channels (channels 1 and 2) for RNA analysis
    %pythonCommand = ['python "I:\OneDrive - California Institute of Technology\Long Cai - 1\script\streamline-seqFISH-master-20210627\streamline-seqFISH-master\src\fiducial_marker_alignment-master\example_matlab_v3rna.py" ' ...
    %    home_dir ' ' ref_fname ' ' ref_fname2 ' ' ro_fname ' ' savefname ' ' posname];
    
    %system(pythonCommand);

    %% get offsets from mat file     
    offsetsName = sprintf(savefnameMat, position);%['20191125_pos' num2str(position) '_offsets_initial_AdjustedPointRatio-logdots-20191213_xyse3_zse3_xyte2_zte2_xyme1_zme1.csv'];
    offsetsPath = getfile(offsetsDir, offsetsName, 'match');
    
    % get points and apply offsets
    offsets = alignpointswrapper_v3lncrna(chArray, [], offsetsPath, ...
        chaTform, numRounds, folderArray, [], usechabboffsets, usephysicaloffsets);

end