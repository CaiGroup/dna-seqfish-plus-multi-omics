% main script to run 100k DNA seqFISH+/sequential immunofluorescence analysis.
% requirement
% Decoded DNA seqFISH++ spots
% initial_fiducial_markers, final_fiducial_markers folders
% segmentation folder with 3D labeles in cellnum.mat in each position folder
% sequential immunofluorescence raw images
% written by Yodai Takei (based on previous versions from the Cai lab)

%% Experiment dependent Variables
addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts');
javaaddpath 'C:\Program Files\MATLAB\R2019a\java\mij.jar';
javaaddpath 'C:\Program Files\MATLAB\R2019a\java\ij.jar';

experimentDir = 'H:\Yodai\100k\2021-04-04-E14-100k-DNAfull-rep2-plate2';
experimentName = '2021-04-04-E14-100k-DNAfull-rep2-IF-final';
experimentLabel = experimentName;

% original DNA seqFISH analysis folder
experimentDir_ref = 'H:\Yodai\100k\2021-03-30-E14-100k-DNAfull-rep2';
experimentLabel_ref = '2021-06-27-E14-100k-rep2-final';
% decoded points folder
save_folder = 'I:\OneDrive - California Institute of Technology\Long Cai - 1\100k\LC1-100k-001-E14-rep2\output';
save_name = 'LC1-100k-001-012-2021-07-06-E14-100k-rep2-output-finalpoints-paints-decoded-onoff-chAll-pos';

posArray = 0:10; % for decoding. needs Pos0 anyway. Don't start from Pos1 or later.
sizeZ = 27; % number of z slices. in the future, this should be automatically calculted.
folderArray = 0:35;

% import chaTform computed from initial fiducial marker image in python
load('H:\Yodai\100k\2021-03-30-E14-100k-DNAfull-rep2\rep1-yodai-rad-all\Matlab_chatform_rad7_ztransxyaff.mat','tform');
chaTformGlobal = tform;
clear tform

%% Other fixed variables for DNA seqFISH+. Some may still need to be changed to compare to.
numChAll = 3;
numRounds = 1; % barcoding rounds
numChannels = length(folderArray); % barcoding pseudo-channels.
superres = 'radial3d'; % 7x7x3 piexels for fitting
hyb1_ref = true; % whether to align hyb1 (true) or initial_fiducial_markers (false)
for i = 1:3
    chaTform{i,1} = {};
end
for i = 1:3
    physicalTform = {};
end
usechabboffsets = false; % use false for the latest version as chromatic shifts are corrected in the very beginning.
usephysicaloffsets = false; % use false.
filtersigma = true;
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
    'numChannels', 'usechabboffsets', 'usephysicaloffsets', 'filtersigma', ...
    'experimentDir_ref', 'experimentLabel_ref');

for position = posArray
    
    %% Step2: Process - grab the raw points from preprocessed images.
    [rawpoints, intensity, pointsCsvName, I] ...
        = seqfishprocess_initialpoints_v2IF(experimentDir, experimentName, experimentLabel, ...
        position, folderArray, superres, chaTformGlobal, filtersigma, numChAll, experimentDir_ref);

    
    %% Step3: Compute tforms for all hybs with fiducial alignment and output aligned/corrected points and offsets
    % define file name here
    
    % for the original folder
    initFolderName = 'initial_fiducial_markers';
    endingString = [initFolderName '-' experimentLabel_ref];
    refInitialSaveName{position+1} = ['ref-points-pos' num2str(position) '-' endingString '.csv'];
    finalFolderName = 'final_fiducial_markers';
    endingString = [finalFolderName '-' experimentLabel_ref];
    refFinalSaveName{position+1} = ['ref-points-pos' num2str(position) '-' endingString '.csv'];
    
    % update line55 in the script for the .py file directory.
    [offsets] = seqfishtformalign_v3IF(experimentDir, experimentLabel, ...
        position, numRounds, numChAll, folderArray, chaTform, ...
        usechabboffsets, usephysicaloffsets, refInitialSaveName{position+1}, pointsCsvName, refFinalSaveName{position+1},...
        experimentDir_ref, experimentLabel_ref, hyb1_ref);
    % pointsch: already barcoded cell format
    % corrpoints: sequential cell format, cell(numFolders, numCh)
    % offsets: reference image/points are hyb1 ch1.
    
    %% Step4 Output images I if necessary.
    Ipaint = seqfish_outputpreprocessedimages_100k_v2IF(experimentDir, experimentName, experimentLabel, I, offsets, chaTform, physicalTform, position, usechabboffsets, usephysicaloffsets, saveprocessedImages);
    % Ipaint is aligned version of I with two channels.
    clearvars I
    
    % Z-score chromosome paint intensity (hyb61-96) per nucleus.
    segPath = fullfile(experimentDir_ref, 'segmentation', ['pos' num2str(position)], 'cellnum.mat');
    load(segPath,'cellnum');
    
    %% output mean raw IF intensity per marker per cell
    getifv2(position, cellnum, Ipaint, numChannels, experimentDir, experimentLabel);
    
    IpaintZscore = ChromPaintIntensities_zscore(cellnum,Ipaint);
    
    clearvars Ipaint
    
    Ifinal = I_per_zslices(IpaintZscore);
    
    clearvars IpaintZscore
    
    %% output z-scored IF intensity per marker per cell
    parfor z = 1:sizeZ
        finalpoints_per_zslices(experimentDir, experimentLabel, position, numChannels, Ifinal{z}, save_folder, save_name, z);
    end
          
end