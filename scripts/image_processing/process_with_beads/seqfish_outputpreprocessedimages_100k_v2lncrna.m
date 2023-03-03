function I = seqfish_outputpreprocessedimages_100k_v2lncrna(experimentDir, experimentName, experimentLabel, I, offsets, tform2ref, tformDapiRot, position, folderArray, saveprocessedImages)
    % output z = 12    

    folderArrayIdx = folderArray + 1;
    numHybs = length(folderArrayIdx);
    numCh = size(I,2);
    saveDirPath = fullfile(experimentDir, 'analysis', experimentLabel, 'processedImages');
    if exist(saveDirPath, 'dir') ~= 7
        mkdir(saveDirPath);
    end
    
    % apply the fiducial tform to the images
    finalOffset = cell(size(I,1), numCh);
    for f = folderArrayIdx
        for ch = 1:numCh
            finalOffset{f,ch} = maketform2(-offsets{ch}.col(f), offsets{ch}.row(f), -offsets{ch}.z(f));
            I{f,ch} = imwarp(I{f,ch}, finalOffset{f,ch}, 'OutputView', imref3d(size(I{f,ch})));
            I{f,ch} = imwarp(I{f,ch}, tform2ref, 'OutputView', imref3d(size(I{f,ch})));
            I{f,ch} = imwarp(I{f,ch}, tformDapiRot, 'OutputView', imref3d(size(I{f,ch})));
        end
    end
    
    %savePath = fullfile(saveDirPath, ['alignedcorr-processed-variables-pos' num2str(position) '-' experimentName '.mat']);
    %save(savePath, 'offsets', 'chaTform', 'physicalTform', 'position', 'usechabboffsets', 'usephysicaloffsets', 'finalOffset');
    
    if saveprocessedImages
        
        Iout = cell(numHybs*3,1);
        Iout(1:numHybs,1) = I(folderArrayIdx,1);
        Iout(numHybs+1:numHybs*2,1) = I(folderArrayIdx,2);
        Iout(numHybs*2+1:numHybs*3,1) = I(folderArrayIdx,3);
        for f = 1:length(Iout)
            Iout{f} = Iout{f}(:,:,12:14); % just for center part of the images for visual check.
        end
        saveTifPath = fullfile(saveDirPath, ['alignedcorr-processed-ch1-ch2-ch3-I-pos' num2str(position) '-' experimentName '.tif']);
        savechannelsimagej(Iout, saveTifPath);

    end
    
    I = I(folderArrayIdx,2:3);
    
    saveImPath = fullfile(saveDirPath, ['alignedcorr-ch2-ch3-I-pos' num2str(position) '.mat']);
    save(saveImPath, 'I', '-v7.3');
    
end

