function I = outputprocessimages_chTcorrectedIFext_v2(experimentDir, ...
    position, folderArray, numCh, typedots, superres, medIntensity, numRefPoints, threshold, refposition, ...
    experimentLabel, roimask, thresholdadjust, filtersigma, chaTform, back, edge)
% output the processeed images depending on the position so the points can
% be extracted for the bead alignment
% chromatic aberration is corrected with the images in the very first step.
    
    %% Declare Variables
    points = cell(length(folderArray),numCh);
    intensity = cell(length(folderArray), numCh);
    intcheck = cell(length(folderArray), numCh);
    sigma = cell(length(folderArray), numCh);
    xyPixSize = 1;
    zPixSize = 1;
    maxXY = 2048;
    maxZ = [];
    min = 1;
    adjustedPoints = cell(length(folderArray),numCh);
    adjustedIntensity = adjustedPoints;
    medianError = adjustedPoints;
    adjustedThreshold = ones(length(folderArray),numCh) * 999999;
    shadingcorr = [];
            


    %% load pos1-4 images, grab points using same treshold vales
    % load threshold
    %load('I:\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH\threshold\thresholdAllCh-pos0-numHybCycles85-numCh-3-E14-DNA-seqFISH+rep3-1-DNAFISH.mat');
    % get images and grab points
    sizeC = [];
    I = cell(length(folderArray), numCh);
    numPoints = I;
    dapiI = cell(length(folderArray), 1);
    allIms = cell(length(folderArray),1);
    folderArrayIdx = folderArray + 1;
    tic
    parfor f = folderArrayIdx
        fprintf('Retrieving Position %.0f Folder %.0f images\n', position, f-1);
        listing = dir([experimentDir '\HybCycle_' num2str(f-1) '\*MMStack_Pos' num2str(position) '.ome.tif']);
        imageName = listing(1).name;
        %imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
        imagePath = fullfile(experimentDir, ['HybCycle_' num2str(f-1)], imageName);
        [allIms{f}, sizeC, sizeZ, ~, ~] = grabimseries(imagePath, position);
    end
    
    % added to correct chromatic aberration in the very beginning.
    for f = folderArrayIdx
        for ch = 1:numCh
            allIms{f}{ch} = imwarp(allIms{f}{ch}, chaTform{ch}, 'OutputView', imref3d(size(allIms{f}{ch})));
        end
    end

    %% Align Dapi for Background Images for All Positions
    %fprintf('Aligning Dapi for Background Images...\n');
    %backgroundFolderName = 'final_fiducial_markers';%'initial_background';
    %backImBasePath = fullfile(experimentDir, backgroundFolderName);
    %listing = dir([experimentDir '\' backgroundFolderName '\*MMStack_Pos' num2str(position) '.ome.tif']);
    %backImPath = fullfile(backImBasePath, listing(1).name);
    %listing = dir([experimentDir '\HybCycle_0\*MMStack_Pos' num2str(position) '.ome.tif']);
    %dapiRefName = listing(1).name;
    %dapiRefPath = fullfile(experimentDir, 'HybCycle_0', dapiRefName);
  
    %shadingcorr = shadingcorrection(backIms(1:numCh));
    for f = folderArrayIdx
        for ch = 1:numCh
            fprintf('folder: %.0f; channel %.0f\n', f, ch);
            % Apply the shading correctionsmean
            %I{f,ch} = uint16(double(allIms{f}{ch}) ./ double(shadingcorr{ch}));
            I{f,ch} = uint16(double(allIms{f}{ch}));

            % ImageJ Rolling Ball Back Subtract to remove noise using rad 3
            % replace with deconvolution - need to test first
            if back
                uniqueString = 'imageTempProcess-90jf03j';
                I{f,ch} = imagejbackgroundsubtraction(I{f,ch}, uniqueString,...
                    experimentDir);
            end
            % find exterior for the IF images
            if edge
                I{f,ch} = imagejfindedges(I{f,ch}, uniqueString, experimentDir);
            end
        end
    end

    toc

end