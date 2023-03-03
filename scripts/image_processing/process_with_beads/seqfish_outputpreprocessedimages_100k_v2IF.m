function I = seqfish_outputpreprocessedimages_100k_v2IF(experimentDir, experimentName, experimentLabel, I, offsets, chaTform, physicalTform, position, usechabboffsets, usephysicaloffsets, saveprocessedImages)
    % output z = 12    

    numHybs = size(I,1);
    numCh = size(I,2);
    saveDirPath = fullfile(experimentDir, 'analysis', experimentLabel, 'processedImages');
    if exist(saveDirPath, 'dir') ~= 7
        mkdir(saveDirPath);
    end
    
    % apply the fiducial tform to the images
    finalOffset = cell(numHybs, numCh);
    for f = 1:numHybs
        for ch = 1:numCh
            finalOffset{f,ch} = maketform2(-offsets{ch}.col(f), offsets{ch}.row(f), -offsets{ch}.z(f));
            I{f,ch} = imwarp(I{f,ch}, finalOffset{f,ch}, 'OutputView', imref3d(size(I{f,ch})));
        end
    end
    
    savePath = fullfile(saveDirPath, ['alignedcorr-processed-variables-pos' num2str(position) '-' experimentName '.mat']);
    save(savePath, 'offsets', 'chaTform', 'physicalTform', 'position', 'usechabboffsets', 'usephysicaloffsets', 'finalOffset');
    
    if saveprocessedImages
        
        Iout = cell(size(I,1)*3,1);
        Iout(1:size(I,1),1) = I(:,1);
        Iout(size(I,1)+1:size(I,1)*2,1) = I(:,2);
        Iout(size(I,1)*2+1:size(I,1)*3,1) = I(:,3);
        for f = 1:length(Iout)
            %Iout{f} = Iout{f}(:,:,12:14); % just for center part of the images for visual check (default for cell culture)
            Iout{f} = Iout{f}(:,:,7:9);
            %Iout{f} = Iout{f}(:,:,26:28); % just for center part of the images for visual check.
        end
        %saveTifPath = fullfile(saveDirPath, ['alignedcorr-processed-ch1-ch2-ch3-I-pos' num2str(position) '-' experimentName '.tif']);
        saveTifPath = fullfile(saveDirPath, ['alignedcorr-processed-ch1-ch2-ch3-I-pos' num2str(position) '-z7-9.tif']);
        savechannelsimagej(Iout, saveTifPath);

    end
    
    I = I(:,1:2); % remove channel 3 for output
        
end

