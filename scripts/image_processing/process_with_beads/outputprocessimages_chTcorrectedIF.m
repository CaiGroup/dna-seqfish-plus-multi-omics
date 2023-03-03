function [points, intensity, shadingcorr, adjustedThreshold, I] = outputprocessimages_chTcorrectedIF(experimentDir, ...
    position, folderArray, numCh, typedots, superres, medIntensity, numRefPoints, threshold, refposition, ...
    experimentLabel, roimask, thresholdadjust, filtersigma, chaTform)
% output the processeed images depending on the position so the points can
% be extracted for the bead alignment
% chromatic aberration is corrected with the images in the very first step.

    saveDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points');
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end

    saveFigDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', ['points_debug_pos' num2str(position)]);
    if exist(saveFigDir, 'dir') ~= 7
        mkdir(saveFigDir);
    end
    
    
    
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
    for f = 1:length(folderArray)
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
    
    
    imageJBackSubtract = false;
    %shadingcorr = shadingcorrection(backIms(1:numCh));
    for f = 1:length(folderArray)
        for ch = 1:numCh
            fprintf('folder: %.0f; channel %.0f\n', f, ch);
            % Apply the shading correctionsmean
            %I{f,ch} = uint16(double(allIms{f}{ch}) ./ double(shadingcorr{ch}));
            I{f,ch} = uint16(double(allIms{f}{ch}));

            if ch == numCh %imageJBackSubtract
                % ImageJ Rolling Ball Back Subtract to remove noise using rad 3
                % replace with deconvolution - need to test first
                uniqueString = 'imageTempProcess-90jf03j';
                I{f,ch} = imagejbackgroundsubtraction(I{f,ch}, uniqueString,...
                    experimentDir);
            end


            if position ~= refposition && thresholdadjust
                %% get threshold if adjusted
                [adjustedThreshold(f,ch), adjustedIntensity{f,ch}, adjustedPoints{f,ch}] = ...
                    outputadjustedthresholdbyone(I{f,ch}, medIntensity{f,ch}, numRefPoints{f,ch}, threshold(f,ch), ...
                    typedots, roimask);
            else
                adjustedThreshold = threshold;
            end



            %% grab the points
            if ch == numCh
                saveFigPath = fullfile(saveFigDir, ['figure-Pos' num2str(position) '-Folder' num2str(f) '-Ch' num2str(ch) '.fig']);
                [pointsTemp, intcheck{f,ch}, ~, ~] = detectdotsv2(I{f, ch}, adjustedThreshold(f,ch), typedots, true, saveFigPath, 1);
                switch superres
                    case 'radial3d'
                        [points{f,ch}, intensity{f,ch}, sigma{f,ch}] = SuperResPoints(pointsTemp,I{f,ch},xyPixSize,zPixSize);
                    case 'gaussian'
                        [points{f,ch}, intensity{f,ch}] = getgaussian(pointsTemp, I{f,ch});
                        [~, ~, sigma{f,ch}] = SuperResPoints(pointsTemp,I{f,ch},xyPixSize,zPixSize, false);
                end

                % filter using sigma
                if filtersigma
                    removeidx = [];
                    for i = 1:size(points{f,ch},1)
                        sigcheck = sigma{f,ch}(i);
                        int = intensity{f,ch}(i);
                        if sigcheck < 0.58 % just to remove shot noise.
                            removeidx = cat(1, removeidx, i);
                        end
                    end
                    points{f,ch}(removeidx, :) = [];
                    intensity{f,ch}(removeidx, :) = [];
                    sigma{f,ch}(removeidx, :) = [];
                end
            end

            %}
            %I{f,ch} = [];

        end
    end

    saveDirPath = fullfile(experimentDir, 'analysis', experimentLabel, 'points');
    if exist(saveDirPath, 'dir') ~= 7
        mkdir(saveDirPath);
    end
    if position ~= refposition
        savePath = fullfile(saveDirPath, ['points-int-thresh-pos' num2str(position) '-' experimentLabel '.mat']);
    else
        savePath = fullfile(saveDirPath, ['points-int-thresh-pos' num2str(position) '-' experimentLabel '.mat']);
    end
    save(savePath, 'points', 'intensity', 'adjustedThreshold', 'shadingcorr', 'sigma');
    toc

end