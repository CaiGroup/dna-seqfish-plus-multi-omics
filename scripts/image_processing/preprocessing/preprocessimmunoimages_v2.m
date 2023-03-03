function preprocessimmunoimages_v2(experimentName, experimentDir, position,  ...
    divideIms, finalFixedPath, finalMovingPath, finalSeqPath, DNAoffsetDir)
% preprocessimages takes raw images from an experiment and preprocesses the
% images: align by dapi, align to reference image, subtract the background, and background subtract
% using the imageJ rolling ball algorithm.
%
% Inputs: experiment name (used for naming files), experiment directory
% (main directory with HybCycle_[i] images, etc.), position or fov, folder array
% (array for number of folders ex. 0:4 for [0,1,2,3,4] folders.
%
% Output: processed images; saves the images in the project folder as a
% tiff file and a mat data file
%
% Requirements: images need at least 4 Z-slices for the alignment
%
% Options: decide to use background images or not - usually used to remove
% tissue background or autofluorescence.
%
% addpath('C:\github\streamline-seqFISH\src\preprocessing\bfmatlab', '-end');
%
% Dependencies: 
% 1. Miji.jar is in the Fiji.app\scripts directory. If not, download and
% install at: ......
% 2. Functions imagejbackgroundsubtraction, grabimseries, shadingcorrection,
% imdivide, grabtform, getdirectory
%
% Time: 3 field of views (without saving hyb processed images) took 5 hours
% - increase speed using parfor loops or increase read/write speed for
% saving.
%
% To Do:
% Add an option to save images; default is on
%
% Options: 
% 1. useBackgroundImages: boolean to use background images for subtraction
% 2. backgroundFolderName: change default folder name in experimentDir
% from 'initial_background' to new folder
% 3. dapiRefPath: path to images that will be used as the dapi reference to
% align images.
% 4. imageJBackSubtract: boolean to use imageJ rolling ball background
% subtraction.
% 5. subtractBackground: boolean to subtract background - useful for high
% bakcground images or autofluorescence.
% 6. saveProcessedHybIms: boolean to save processed Hyb images
% default Values for optional arguments:
% 7. divideIms: divides dapi images into 4 to use gradient descent
% for 3d alignment - good for stacks with less than 16 zslices
% [true, 'initial_background', [], false, false, false]; 
%
% added rotation correction and downsampling option by YT.
%
%
% 
% Date: 07/09/2021
% Author: Yodai Takei (modified original script from Nico Pierson)
% ytakei@caltech.edu

    tic
    %% Check if Fiji and bfmatlab is in the path
    fijiDirectory = checkfijipath();
    checkbfmatlabpath();
    
    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    endingDateString = datestr(dateStart, formatDate);
    
    
    %% Initialize Variables
    folderSaveName = 'processedimages_rotalign';


    % Make save directory: default 'organizehybs' folder
    saveDir = fullfile(experimentDir, folderSaveName, ['pos' num2str(position)]);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    
    
    %% Get the tform for RNA to DNA experiment
    imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
    imFixedPath = fullfile(finalFixedPath, imageName);
    imMovingPath = fullfile(finalMovingPath, imageName);
    imSeqPath = fullfile(finalSeqPath, imageName);
    % initial RNA fiducial marker image without DAPI (3 channels)
    [DNAIms, sizeCDNA, sizeZDNA, ~, ~] = grabimseries(imFixedPath, position);
    % final alignment image including DAPI (4 channels)
    [RNAIms, sizeCRNA, sizeZRNA, ~, ~] = grabimseries(imMovingPath, position);
    % DNA seqFISH+ image
    [SeqIms, sizeCSeq, sizeZSeq, ~, ~] = grabimseries(imSeqPath, position);
    if ~isempty(DNAoffsetDir)
        listing2 = dir(fullfile(DNAoffsetDir,['pos' num2str(position)],'imagesHybDapi-*.mat'));
        load(fullfile(listing2(1).folder,listing2(1).name),'tformDapi');
        DNADapi_Im = imwarp(SeqIms{sizeCSeq}, tformDapi{1,1}, 'OutputView', imref3d(size(SeqIms{sizeCSeq})));
    else
        DNADapi_Im = SeqIms{sizeCSeq};
    end
    
    initialRadius = 0.0625; %0.0625 for 3d is default
    numIterations = 100; % 100 is default
    
    % compute tform between the 100k RNA initial fiducial and final segmentation
    % alignment with channel 3.
    if divideIms
        [RNAIms4, ~] = imdivideby4(RNAIms{sizeCDNA});
        [DNAIms4, ~] = imdivideby4(DNAIms{sizeCDNA});
        % multimodal alignment instead of monomodal alignment
        tform2ref_pre = grabtformRot(RNAIms4, DNAIms4,initialRadius,numIterations);
    else
        tform2ref_pre = grabtformRot(RNAIms{sizeCDNA}, DNAIms{sizeCDNA},initialRadius,numIterations);
    end
    
    % remove rotation fraction for these images.
    tform2ref_pre.T(1,2) = 0; 
    tform2ref_pre.T(2,1) = 0; 
    
    % align final alignment DAPI to initial fiducial marker image
    % then compute the translation and rotation tform between the RNA initial fiducial image and DNA hyb1 image 
    
    if divideIms
        % align rna final segmentation alignment image dapi to rna initial fiducial marker image.
        newIm_align_fid = imwarp(RNAIms{sizeCDNA}, tform2ref_pre, 'OutputView', imref3d(size(RNAIms{sizeCDNA})));
        newIm_align_pre = imwarp(RNAIms{sizeCRNA}, tform2ref_pre, 'OutputView', imref3d(size(RNAIms{sizeCRNA})));
        % compute and align rna initial fiducial marker image to hyb1 DNA image
        [newIm_align_pre4, ~] = imdivideby4(newIm_align_pre);
        [DNADapi_Im4, ~] = imdivideby4(DNADapi_Im);
        tform2ref = grabtform(newIm_align_pre4, DNADapi_Im4,initialRadius,numIterations);
        newIm_align = imwarp(newIm_align_pre, tform2ref, 'OutputView', imref3d(size(newIm_align_pre)));
        [newIm_align4, ~] = imdivideby4(newIm_align);
        tformDapiRot = grabtformRot(newIm_align4, DNADapi_Im4, initialRadius, numIterations);
        tformDapiRot.T(4,3) = 0; 
        newIm_align_final = imwarp(newIm_align, tformDapiRot, 'OutputView', imref3d(size(newIm_align)));
    else
        % align rna final segmentation alignment image dapi to rna initial fiducial marker image.
        newIm_align_fid = imwarp(RNAIms{sizeCDNA}, tform2ref_pre, 'OutputView', imref3d(size(RNAIms{sizeCDNA})));
        newIm_align_pre = imwarp(RNAIms{sizeCRNA}, tform2ref_pre, 'OutputView', imref3d(size(RNAIms{sizeCRNA})));
        % compute and align rna initial fiducial marker image to hyb1 DNA image
        tform2ref = grabtform(newIm_align_pre, DNADapi_Im,initialRadius,numIterations);
        newIm_align = imwarp(newIm_align_pre, tform2ref, 'OutputView', imref3d(size(newIm_align_pre)));
        tformDapiRot = grabtformRot(newIm_align, DNADapi_Im, initialRadius, numIterations);
        tformDapiRot.T(4,3) = 0; 
        newIm_align_final = imwarp(newIm_align, tformDapiRot, 'OutputView', imref3d(size(newIm_align)));
    end
    
    I = cell(5,1);
    I{1,1} = DNAIms{sizeCDNA}; % RNA FISH initial fiducial marker
    I{2,1} = newIm_align_fid; % RNA FISH final segmentation fiducial
    I{3,1} = newIm_align; % RNA FISH DAPI before rotation
    I{4,1} = newIm_align_final; % RNA FISH DAPI after rotation
    I{5,1} = DNADapi_Im; % DNA FISH DAPI
    
    savePath = fullfile(saveDir, ['preProcessedData-pos' num2str(position) '-' experimentName '-' endingDateString '.mat']);
    save(savePath, 'position', 'tform2ref','tformDapiRot');
    
    savePath2 = fullfile(saveDir,'alinged_rotated_DAPI.tif');
    savechannelsimagej(I, savePath2);
    
    toc
end