function [] = decodeallcellshpc(experimentDir, experimentName, position, ...
        numRounds, numChannels, points, cell, saveDir, iter, offpercent, varargin)
% decodes all the cells from the input of points provided and
% finding the false positive rate
% 
% Requirements: 
% 1. points
% 2. thresholds (round x channel matrix) - in 'threshold' folder in project
% directory
% 3. barcode key in excel - in 'barcodekey' folder in project directory
% 4. RoiSet.zip for 2d segmentation - in project directory
%
% Inputs: project directory path, project name will be used to save data
% and directories, position (aka field of view), number of barcode rounds,
% number of channels, points, cell number, save Directory to save HPC-Data,
% number of iteration (usually there is more than 1 iteration.
%
% Requirements: 'segmentation' folder for rois, points
%
% Dependencies: process package
%
% Optional Parameters: 4. allowed diff 5. square
% root radius
%
% Default: 1, 6
%
% To do: 
% 1. Add option to find points - how if the images are so large
% 2. Add segmentation for labeled images from ilastik
% 3. 3d functions for segmentation
% 5. Update Readme.txt and make it clear and concise
% 6. Add an example
%
% Date: 9/9/2019
% Author: Nico Pierson

    %% Set up optional Parameters
    argsLimit = 7;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:processimages:TooManyInputs', ...
            'requires at most 9 optional inputs');
    end
    % Error for type of arguments
    if numvarargs > 0
        if ~isnumeric(varargin{1}) && ~isscalar(varargin{1})
            error('src:processimages:WrongInput', ...
                'processmages sqrtradius requires type int');
        end
    end
    if numvarargs > 1
        if ~isnumeric(varargin{2}) && ~isscalar(varargin{2})
            error('src:processimages:WrongInput', ...
                'processimages var alloweddiff requires type int');
        end
    end
    if numvarargs > 2
        if ~isnumeric(varargin{3}) && ~isscalar(varargin{3})
            error('src:processimages:WrongInput', ...
                'processimages var channel requires type int');
        end
    end
    if numvarargs > 3
        if ~isnmatrix(varargin{4})
            error('src:processimages:WrongInput', ...
                'processimages var minseeds requires type mat');
        end
    end
    if numvarargs > 4
        if varargin{5} ~= 0 && varargin{5} ~= 1
            error('src:processimages:WrongInput', ...
                'processimages var removePoints requires type matrix');
        end
    end
    if numvarargs > 5
        if isempty(varargin{6})
            error('src:processimages:WrongInput', ...
                'processimages var dotlocations is required');
        end
    end
    if numvarargs >= argsLimit
        if ~ismatrix(varargin{7})
            error('src:processimages:WrongInput', ...
                'processimages var removeInd requires type matrix');
        end
    end
    % set defaults for optional inputs
    optargs = {6, 1, [], [], false, [], []};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [sqrtradius, alloweddiff, channel, minseeds, removePoints, dotlocations, removeInd] = optargs{:}; % beadPoints used to remove bead points

    
    
    %% Read the Barcode excel file
    barcodeFolder = 'barcodekey';
    if isempty(channel)
        barcodekeyPath = getfile(fullfile(experimentDir,barcodeFolder), 'barcodekey');
    else
        channelFolder = ['ch' num2str(channel)];
        barcodekeyPath = getfile(fullfile(experimentDir,barcodeFolder, channelFolder), 'barcodekey');
    end
    barcodekey = readbarcode(barcodekeyPath, 'no header');
    % remove 1st of numbers from 1st column in Linus's data
    barcodekey.barcode(:,1) = [];
    % generate all possible barcodes
    offtarget = generatebarcodekeyseqfish();
    % remove the on-target barcodes from the off-target list
    for i = 1:size(barcodekey.barcode, 1)
        % find the index of the row
        [tf, index] = ismember(barcodekey.barcode(i,:), offtarget.barcode, 'rows');

        % delete the point
        if tf
            offtarget.barcode(index,:) = [];
            offtarget.names(index) = [];
        end
    end
    offtargetString = 'offtarget_';
    offtarget.names = strcat(offtargetString, offtarget.names);
    numOffTargetBarcodes = size(offtarget.barcode, 1);


 

    
    if removePoints && ~isempty(dotlocations)
        %% Remove Points - make function to remove points within a certain radius
        % need the previous points and the dotlocations, with indices to remove
        % certain points
        points = removepoints(removePoints, dotlocations, removeInd);
    end
    
    
    
    %% Decode Points for each ROI or labeled cell
    % Get the path for segmentation - need to make useful for 3d
    segmentFolder = 'segmentation';
    segmentPath = fullfile(experimentDir, ['RoiSet_Pos' num2str(position) '.zip']);
    %segmentPath = fullfile(experimentDir, segmentFolder, ['Pos' num2str(position)], 'RoiSet.zip');
    segment = 'roi'; % can use ['roi', '3d']
    savePath = fullfile(saveDir, 'hpc-data');
    if exist(savePath, 'dir') ~= 7
        mkdir(savePath);
    end
    
    
    
    %% Add 10% of off-target barcodes to end of list
    % May need to loop over a few times for this type of check
    offtargetDivideFactor = offpercent / 100;
    numOffTargets = ceil(numOffTargetBarcodes * offtargetDivideFactor);
    rng('shuffle');
    idxOff = sort(randperm(numOffTargetBarcodes, numOffTargets));
    % add to the end of the list for on target barcodes
    barcodekeyfull.names = cat(1, barcodekey.names, offtarget.names(idxOff));
    barcodekeyfull.barcode = cat(1, barcodekey.barcode, offtarget.barcode(idxOff,:));
    
    
    
    %% decode each cell in a position and save it
    decodecellhpc(points, experimentName, position, channel, sqrtradius, ...
    alloweddiff, numRounds, numChannels, cell, savePath, segmentPath, segment, ...
    barcodekeyfull, iter, offpercent, minseeds)

end