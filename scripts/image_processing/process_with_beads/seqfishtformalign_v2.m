function [pointsch, offsets, corrpoints] = seqfishtformalign_v2(experimentDir, experimentLabel, ...
    position, numRounds, numCh, ...
    folderArray, chaTform, physicalTform, usechabboffsets, usephysicaloffsets, ...
    refSaveName, hybSaveName, refSaveName2)   
    
    % 2nd phase for decoding points
    %addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
    %addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts\', '-end');
    %addpath('C:\github\streamline-seqFISH\src\beadalignment\', '-end');


    
    %% variables
    chArray = 1:numCh;
    numHybs = length(folderArray);
    pointsDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points');%, 'pre_formated');
    preformatDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'pre_formated');
    if exist(preformatDir, 'dir') ~= 7
        mkdir(preformatDir);
    end
    offsetsDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'positions');
    if exist(offsetsDir, 'dir') ~= 7
        mkdir(offsetsDir);
    end
    
    pos = position;
    hybSaveNameNew = ['hybridization-points-pos' num2str(position) '.csv'];
    hybCsvOldPath = fullfile(preformatDir, hybSaveName);
    hybCsvNewPath = fullfile(preformatDir, sprintf(hybSaveNameNew, pos));
    copyfile(hybCsvOldPath, hybCsvNewPath);
    refSaveNameNew = ['reference-points-pos' num2str(position) '.csv'];
    refCsvOldPath = fullfile(preformatDir, refSaveName);
    refCsvNewPath = fullfile(preformatDir, sprintf(refSaveNameNew, pos));
    copyfile(refCsvOldPath, refCsvNewPath);
    
    refSaveNameNew2 = ['reference2-points-pos' num2str(position) '.csv'];
    refCsvOldPath2 = fullfile(preformatDir, refSaveName2);
    refCsvNewPath2 = fullfile(preformatDir, sprintf(refSaveNameNew2, pos));
    copyfile(refCsvOldPath2, refCsvNewPath2);



    %% Align using Beads
    home_dir = strcat('"', pointsDir, '"');
    ref_fname = strcat('"', fullfile(preformatDir,refSaveNameNew), '"');
    ref_fname2 = strcat('"', fullfile(preformatDir,refSaveNameNew2), '"');
    ro_fname = strcat('"', fullfile(preformatDir,hybSaveNameNew) , '"');
    savefnameMat = ['point-offsets_pos' num2str(position) '_'];
    savefname = strcat('"', sprintf(savefnameMat, position), '"');
    posname = strcat('"', num2str(position), '"');

    % Run the python comand for bead alignment
    
    % new version for more robust alignment correction.
    pythonCommand = ['python "I:\OneDrive - California Institute of Technology\Long Cai - 1\script\streamline-seqFISH-master-20210627\streamline-seqFISH-master\src\fiducial_marker_alignment-master\example_matlab_v3.py" ' ...
        home_dir ' ' ref_fname ' ' ref_fname2 ' ' ro_fname ' ' savefname ' ' posname];
    
    system(pythonCommand);

    %% get offsets from mat file     
    offsetsName = sprintf(savefnameMat, position);%['20191125_pos' num2str(position) '_offsets_initial_AdjustedPointRatio-logdots-20191213_xyse3_zse3_xyte2_zte2_xyme1_zme1.csv'];
    pointsName = sprintf('points-int-thresh-pos%d', position); % need to grab initial points
    offsetsPath = getfile(offsetsDir, offsetsName, 'match');
    pointsPath = getfile(pointsDir, pointsName, 'match');

    % run alignment code until offsets is output
    limit = 10;
    limit_iter = 1;
    while(exist(offsetsPath, 'file') ~= 2)
        % get offsets
        system(pythonCommand);
        offsetsPath = getfile(offsetsDir, offsetsName, 'match');
        limit_iter = limit_iter + 1;
        if limit_iter >= limit
            error 'alignmnet failed to output offsets';
        end
    end
    
    % get points and apply offsets
    [pointsch, offsets] = alignpointswrapper_v3(chArray, pointsPath, offsetsPath, ...
        chaTform, numRounds, folderArray, physicalTform{position+1,1}, usechabboffsets, usephysicaloffsets); % updated to v3
    
    numChannels = round(numHybs / numRounds);% barcoding pseudo-channels.
    corrpoints = cell(numHybs,numCh);
    for ch = chArray % actual imaging channels.
        for ro = 1:numRounds
            for bch = 1:numChannels
                corrpoints{(ro-1)*numChannels+bch,ch} = pointsch{ch,1}{ro,1}(bch).channels;
            end
        end
    end
    

end