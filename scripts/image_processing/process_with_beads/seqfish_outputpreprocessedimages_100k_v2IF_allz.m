function seqfish_outputpreprocessedimages_100k_v2IF_allz(experimentDir, experimentName, experimentLabel, I, offsets, chaTform, physicalTform, position, usechabboffsets, usephysicaloffsets, saveprocessedImages)
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
        
        Iout = cell(size(I,1)*2,1);
        Iout(1:size(I,1),1) = I(:,1);
        Iout(size(I,1)+1:size(I,1)*2,1) = I(:,2);
        
        % selected 8 images
        Iout2 = cell(9,1);
        Iout2{1} = Iout{38}; % SF3A66
        Iout2{2} = Iout{7}; % RNAPIISer5-P
        Iout2{3} = Iout{16}; % H3K27ac
        Iout2{4} = Iout{17}; % H3K36ac
        Iout2{5} = Iout{40}; % H3K4me2
        Iout2{6} = Iout{52}; % MajSat
        Iout2{7} = Iout{19}; % H4K20me3
        Iout2{8} = Iout{46}; % H3K27me3
        Iout2{9} = Iout{33}; % LaminB1
        
        saveTifPath = fullfile(saveDirPath, ['alignedcorr-processed-ch1-ch2-I-pos' num2str(position) '-' experimentName '-all-z.tif']);
        %saveTifPath = fullfile(saveDirPath, ['alignedcorr-processed-ch1-ch2-ch3-I-pos' num2str(position) '-z9-11.tif']);
        savechannelsimagej(Iout2, saveTifPath);

    end
        
end

