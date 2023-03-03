function I = seqfish_outputpreprocessedimages_all(experimentDir, experimentName, experimentLabel, I, offsets, chaTform, physicalTform, position, usechabboffsets, usephysicaloffsets, saveprocessedImages)
    %%%
    % output chromatic aberration corrected "I"
    % save "I", instead of I1, I2, I3.
    % save I1, I2, I3 for z = 10:15 as tif images to visual check.
    % Last edit on 02/13/20 by Yodai Takei
    %%%

    numHybs = size(I,1);
    numCh = size(I,2);
    saveDirPath = fullfile(experimentDir, 'analysis', experimentLabel, 'processedImages');
    if exist(saveDirPath, 'dir') ~= 7
        mkdir(saveDirPath);
    end
    
    % get the tform for channel 1, 2 and 3
    if length(chaTform) < numCh
        numCh = length(chaTform);
    end
    chabbAndPhysicalTform = cell(1, numCh);
    if usechabboffsets && usephysicaloffsets
        for ch = 1:numCh
            chabbAndPhysicalTform{ch} = addtform(chaTform{ch}, physicalTform{ch});
        end
    elseif usechabboffsets
        for ch = 1:numCh
            chabbAndPhysicalTform{ch} = chaTform{ch};
        end
    else
        error('check chromatic aberration');
    end
    
    % apply the final tform to the images
    fiduciaryOffsets = cell(numHybs, numCh);
    finalOffset = cell(numHybs, numCh);
    for f = 1:numHybs
        for ch = 1:numCh
            fiduciaryOffsets{f,ch} = maketform2(-offsets{ch}.col(f), offsets{ch}.row(f), -offsets{ch}.z(f));
            finalOffset{f,ch} = addtform(chabbAndPhysicalTform{ch}, fiduciaryOffsets{f,ch});
            I{f,ch} = imwarp(I{f,ch}, finalOffset{f,ch}, 'OutputView', imref3d(size(I{f,ch})));
        end
    end
    
    if saveprocessedImages
        % save the processed and aligned images for ch1, ch2 (and ch3).
        saveImName = ['alignedcorr-processed-I-pos' num2str(position) '-' experimentName '.mat'];
        saveImPath = fullfile(saveDirPath, saveImName);
        save(saveImPath, 'I', '-v7.3');
        
        I1 = I(:,1);
        for f = 1:length(I1)
            I1{f} = I1{f}(:,:,10:15); % just for center part of the images for visual check.
        end
        saveTifPath = fullfile(saveDirPath, ['alignedcorr-processed-ch1-I-pos' num2str(position) '-' experimentName '.tif']);
        savechannelsimagej(I1, saveTifPath);
        clearvars I1
        
        
        I2 = I(:,2);
        for f = 1:length(I2)
            I2{f} = I2{f}(:,:,10:15); % just for center part of the images for visual check.
        end
        saveTifPath = fullfile(saveDirPath, ['alignedcorr-processed-ch2-I-pos' num2str(position) '-' experimentName '.tif']);
        savechannelsimagej(I2, saveTifPath);
        clearvars I2
        
        if numCh > 2
            I3 = I(:,3);
            for f = 1:length(I3)
                I3{f} = I3{f}(:,:,10:15); % just for center part of the images for visual check.
            end
            saveTifPath = fullfile(saveDirPath,['alignedcorr-processed-ch3-I-pos' num2str(position) '-' experimentName '.tif']);
            savechannelsimagej(I3, saveTifPath);
            clearvars I3
        end
    
        savePath = fullfile(saveDirPath, ['alignedcorr-processed-variables-pos' num2str(position) '-' experimentName '.mat']);
        save(savePath, 'offsets', 'chaTform', 'physicalTform', 'position', 'usechabboffsets', 'usephysicaloffsets', 'finalOffset');
    end

end

