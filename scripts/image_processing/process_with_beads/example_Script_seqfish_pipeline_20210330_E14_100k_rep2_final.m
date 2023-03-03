% main script to run 100k DNA seqFISH++ analysis.
% requirement
% HybCycle_0 to HybCycle_95 folders for 100k DNA seqFISH++ images
% initial_fiducial_markers, final_fiducial_markers folders
% segmentation folder with 3D labeles in cellnum.mat in each position folder
% thretholding values in threshold folder as csv file
% written by Yodai Takei (based on previous versions from the Cai lab)

%% Experiment dependent Variables
addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts');
javaaddpath 'C:\Program Files\MATLAB\R2019a\java\mij.jar';
javaaddpath 'C:\Program Files\MATLAB\R2019a\java\ij.jar';

experimentDir = 'H:\Yodai\100k\2021-03-30-E14-100k-DNAfull-rep2';
experimentName = '2021-06-27-E14-100k-rep2-final';
experimentLabel = experimentName;
posArray = 0:10; % for decoding. needs Pos0 anyway. Don't start from Pos1 or later.
posArrayAll = 0:10;% for fiducial markers.
sizeZ = 27; % number of z slices. in the future, this should be automatically calculted.

% import chaTform computed from initial fiducial marker image in python
load('H:\Yodai\100k\2021-03-30-E14-100k-DNAfull-rep2\rep1-yodai-rad-all\Matlab_chatform_rad7_ztransxyaff.mat','tform');
chaTformGlobal = tform;
clear tform

%% Other fixed variables for DNA seqFISH+. Some may still need to be changed to compare to.
numChAll = 3;
folderArray = 0:95;
numRounds = 1; % barcoding rounds
numChannels = 96; % barcoding pseudo-channels.
numpointshyb = 60; % loci imaging hybs
superres = 'radial3d'; % 7x7x3 piexels for fitting
chaTform = {};
chaTform{1,1} = {};
chaTform{2,1} = {};
chaTform{3,1} = {};
usechabboffsets = false; % use false for the latest version as chromatic shifts are corrected in the very beginning.
usephysicaloffsets = false; % use false.
filtersigma = true;
segment = 'cellnum'; % this should be cellnum with current version.
saveprocessedImages = true;


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
% chromatically correct images first and then save the points
[~, ~, physicalTform, refPointsAlignedInitial, refInitialSaveName] ...
    = debugbeads_cellnum_v2(experimentDir, posArrayAll, numChAll, sizeZ, superres, experimentLabel,chaTformGlobal);

for position = posArray
    
    %% Step2: Process - grab the raw points from preprocessed images.
    [rawpoints, intensity, pointsCsvName, I] ...
        = seqfishprocess_initialpoints_v2(experimentDir, experimentName, experimentLabel, ...
        position, folderArray, superres, chaTformGlobal, filtersigma, numChAll);

    
    %% Step3: Compute tforms for all hybs with fiducial alignment and output aligned/corrected points and offsets
    % define file name here
    refFinalSaveName{position+1} = strrep(refInitialSaveName{position+1},'initial_fiducial','final_fiducial');
    % update line55 in the script for the .py file directory.
    [pointsch, offsets, corrpoints] = seqfishtformalign_v3(experimentDir, experimentLabel, ...
        position, numRounds, numChAll, folderArray, chaTform, physicalTform, ...
        usechabboffsets, usephysicaloffsets, refInitialSaveName{position+1}, pointsCsvName, refFinalSaveName{position+1});
    % pointsch: already barcoded cell format
    % corrpoints: sequential cell format, cell(numFolders, numCh)
    % offsets: reference image/points are hyb1 ch1.
    
    %% Step4 Output images I if necessary.
    Ipaint = seqfish_outputpreprocessedimages_100k_v2(experimentDir, experimentName, experimentLabel, I, offsets, chaTform, physicalTform{position+1}, position, usechabboffsets, usephysicaloffsets, saveprocessedImages);
    
    %% Step5: decoding for the ch1-3 barcoding dataset.
    clearvars I
    
    % Z-score chromosome paint intensity (hyb61-96) per nucleus.
    segPath = fullfile(experimentDir, 'segmentation', ['pos' num2str(position)], 'cellnum.mat');
    load(segPath,'cellnum');
    IpaintZscore = ChromPaintIntensities_zscore(cellnum,Ipaint);
    
    % formatting FISH spots (hyb1-60).
    finalpoints = seqfishformatforch123(pointsch, numpointshyb);
    
    % save all the variables for decoding to change the decoding paramters
    % at the next run if necessary.
    saveDecodeDir = fullfile(experimentDir, 'analysis', experimentLabel, 'decoding_variables');
    if exist(saveDecodeDir, 'dir') ~= 7
        mkdir(saveDecodeDir);
    end
    saveDecodePath = fullfile(saveDecodeDir,['decoding-variables-pos' num2str(position) '-' experimentName '.mat']);
    save(saveDecodePath, 'experimentDir', 'experimentName', 'experimentLabel', 'superres', ...
    'posArray', 'position','folderArray', 'numChAll', 'sizeZ', 'numRounds', ...
    'numChannels', 'usechabboffsets', 'usephysicaloffsets', ...
    'filtersigma', 'pointsch','cellnum', ...
    'finalpoints', 'segment', 'chaTform','offsets','Ipaint','IpaintZscore','-v7.3');

    % formatting to get paint intensity per hyb1-60 points.
    % output as csv file.
    % final decoding should be done with a separate jupyter notebook.
    finalpoints_paints = ChromPaintIntensities_perpoints_v5(experimentDir, experimentName, experimentLabel, finalpoints, IpaintZscore, position, segment);
    output_finalpoints_paints(experimentDir,experimentLabel,experimentName,position,finalpoints_paints);
    
    clearvars Ipaint
    clearvars IpaintZscore

end