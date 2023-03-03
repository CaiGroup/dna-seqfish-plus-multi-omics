function I = seqfishprocess_initialpoints_v2IFext(experimentDir, experimentName, experimentLabel, ...
    position, folderArray, superres, chaTform, filtersigma, numCh, experimentDir_ref)

    % Chromatic aberraction correction was performed on the images in the
    % very first step after extracting the images.
    
    %% Paths
    addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
    addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts\', '-end');
    

    %% Variables
    refposition = 0;
    typedots = 'log';%'exons';% 'log'; % new log filter
    %numCh = 2;
    thresholdadjust = false;


    %% get the reference threshold
    chArray = 1:numCh;
    % load Pos0 most recent thresholds
    %csvPath = fullfile(experimentDir, 'threshold', '002-016-E14-rep2-2-DNAFISH-ch1-2-hyb1-80-threshold-logfilter-Pos0vs3-20191208.csv');
    csvPath = getfile(fullfile(experimentDir_ref, 'threshold'), 'threshold', 'match');
    %threshold = threshcsv2matyodai(csvPath, chArray);
    fileExt = csvPath(end-2:end);
    if strcmp(fileExt, 'mat')
        % just load the data
        t = load(csvPath);
        threshold = t.threshold;
    elseif strcmp(fileExt, 'csv')
        threshold = threshcsv2matyodai(csvPath, chArray);
    end

    %% Output points form the images
    numPointChannels = size(threshold,2);
    I = outputprocessimages_chTcorrectedIFext(experimentDir, ...
        position, folderArray, numCh, typedots, superres, [], [], threshold, ...
        refposition, experimentLabel, [], thresholdadjust, filtersigma, chaTform);

end

