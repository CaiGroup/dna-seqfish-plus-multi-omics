% main script to run RNA seqFISH+ analysis.
% Channel assignmets are optimized for the DNA seqFISH++ experiment scheme.
% written by Yodai Takei

%% Experiment dependent Variables
addpath('C:\Users\Long Cai - 2\Documents\Fiji.app\scripts');
javaaddpath 'C:\Program Files\MATLAB\R2019a\java\mij.jar';
javaaddpath 'C:\Program Files\MATLAB\R2019a\java\ij.jar';

experimentDir = 'G:\Yodai\100k\2021-05-05-cerebellum-RNAfull-rep1';
experimentName = '2021-05-05-cerebellum-RNAfull-rep1-final';
experimentLabel = experimentName;
posArray = 0:2; % for decoding. needs Pos0 anyway. Don't start from Pos1 or later.
posArrayAll = 0:2;% for fiducial markers.
sizeZ = 49; % number of z slices. in the future, this should be automatically calculted.

% import chaTform computed from tetraspeck beads
load('G:\Yodai\100k\Matlab_chatform_rad7_ztransxyaff.mat','tform');
chaTform = tform;
clear tform

%% Other fixed variables for DNA seqFISH+. Some may still need to be changed to compare to.
numChAll = 2;
channels = 1:numChAll;
folderArray = 0:59;
numRounds = 1; % barcoding rounds
numChannels = 60; % barcoding pseudo-channels.
numpointshyb = 60; % loci imaging hybs
superres = 'radial3d'; % 'radial3d', 'gaussian'
usechabboffsets = true; % should be always true for RNA analysis.
usephysicaloffsets = false; % use false.
filtersigma = true;
segment = 'cellnum'; % this should be cellnum with current version.
saveprocessedImages = true;

% for intron seqFISH+
sqrtradius = 3;
alloweddiff = 1;
minseeds = 4;

%% save global variables
saveGlobalDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'variables');
if exist(saveGlobalDir, 'dir') ~= 7
    mkdir(saveGlobalDir);
end
dateStart = datetime;
formatDate = 'yyyy-mm-dd';
dateString = datestr(dateStart, formatDate);
saveGlobalName = ['global-variables-seqfish-pipeline_' experimentLabel '_' dateString];
saveGlobalPath = fullfile(saveGlobalDir, saveGlobalName);
save(saveGlobalPath, 'experimentDir', 'experimentName', 'experimentLabel', 'superres', ...
    'posArray', 'folderArray', 'numChAll', 'sizeZ', 'numRounds', ...
    'numChannels', 'usechabboffsets', 'usephysicaloffsets', 'filtersigma', 'segment');

%% Step1-6: for DNA seqFISH+ analysis.
%% Step1: output fiducial marker points for all positions
% manually threshold the fiducial marker images.
[~, ~, physicalTform, refPointsAlignedInitial, refInitialSaveName] ...
    = debugbeads_cellnum(experimentDir, posArrayAll, numChAll, sizeZ, superres, experimentLabel);

for position = posArray
    
    %% Step2: Process - grab the raw points from preprocessed images.
    [rawpoints, intensity, pointsCsvName, I] ...
        = seqfishprocess_initialpoints(experimentDir, experimentName, experimentLabel, ...
        position, folderArray, superres, chaTform, filtersigma, numChAll);

    %% Step3: Compute tforms for all hybs with fiducial alignment and output aligned/corrected points and offsets
    % define file name here
    refFinalSaveName{position+1} = strrep(refInitialSaveName{position+1},'initial_fiducial','final_fiducial');
    % update line55 in the script for the .py file directory.
    numRounds = 1;
    [pointsch, offsets, corrpoints] = seqfishtformalign_v2rna(experimentDir, experimentLabel, ...
        position, numRounds, numChAll, folderArray, chaTform, physicalTform, ...
        usechabboffsets, usephysicaloffsets, refInitialSaveName{position+1}, pointsCsvName, refFinalSaveName{position+1});
    % pointsch: already barcoded cell format
    % corrpoints: sequential cell format, cell(numFolders, numCh)
    % offsets: reference image/points are hyb1 ch1.
    
    rotDir = dir(fullfile(experimentDir,'processedimages_rotalign',['pos' num2str(position)],'*.mat'));
    load(fullfile(rotDir(1).folder,rotDir(1).name),'tform2ref','tformDapiRot');
    
    for channel = 1:numChAll
        for idx = 1:length(folderArray)
            pointsch{channel,1}{1,1}(idx).channels = transformPointsForward(tform2ref,  pointsch{channel,1}{1,1}(idx).channels);
            pointsch{channel,1}{1,1}(idx).channels = transformPointsForward(tformDapiRot,  pointsch{channel,1}{1,1}(idx).channels);
        end
    end
    clearvars I
    
    % load segmentation label
    segPath = fullfile(experimentDir, 'segmentation', ['pos' num2str(position)], 'cellnum.mat');
    load(segPath,'cellnum');
    
    % non-barcoded mRNA seqFISH analysis
    channel = 1;
    numRounds = 1;
    numChannels = 60;
    points = seqfishformatforchrna(pointsch, channel, numRounds, numChannels, folderArray(1:numRounds*numChannels));
    
    %sqrtradius = 1;% won't be used
    %alloweddiff = 1;% won't be used
    %minseeds = 1;% won't be used
    processimagespoints_sequential(experimentDir, experimentName, ...
    position, numRounds, numChannels, points, segment, sqrtradius, alloweddiff, ...
    channel, minseeds, experimentLabel);

    % in case for dots without segmentaiton for visualization purposes
    processimagespoints_sequential(experimentDir, experimentName, ...
    position, numRounds, numChannels, points, 'whole', sqrtradius, alloweddiff, ...
    channel, minseeds, [experimentLabel '-whole']);

    % intron seqFISH+ analysis
    channel = 2;
    numRounds = 5;
    numChannels = 12;
    points = seqfishformatforchrna(pointsch, channel, numRounds, numChannels, folderArray(1:numRounds*numChannels));
    
    processimagespoints(experimentDir, experimentName, ...
    position, numRounds, numChannels, points, segment, sqrtradius, alloweddiff, ...
    channel, minseeds, experimentLabel);

end