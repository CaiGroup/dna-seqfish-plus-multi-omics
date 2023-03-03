function getifv2(position, cellnum, Ipaint, numChannels, experimentDir, experimentLabel)
% function to output mean and quantile IF intensity per marker per cell as a csv file
% cellnum: labeled image
% Ipaint: aligned IF images
% numChannels: total hyb number
% Immunokey file should be in "immunokey" folder and sorted by hyb and channel
% and the first column should be antibody names.
% Date: 12/17/2021
% Author: Yodai Takei
    
    saveGlobalDir = fullfile(experimentDir,'analysis',experimentLabel,'IF-intensity-per-cell');
    if exist(saveGlobalDir, 'dir') ~= 7
        mkdir(saveGlobalDir);
    end
    
    % load sequential codebook
    immunoPath = getfile(fullfile(experimentDir, 'immunokey'), '.', 'match');
    if exist(immunoPath, 'file') ~= 2
        immunokey = [];
    else
        immunokey = readsequential(immunoPath, 'header');
    end
    
    listSavePath = fullfile(saveGlobalDir, ['output-IF-intensity-per-cell-pos' num2str(position) '.csv']);
    fileID = fopen(listSavePath,'w');
    fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
                'fov', 'hyb','channel','target','cellID', 'mean_intensity',...
                '25th_intensity','median_intensity','75th_intensity');
            
    s = regionprops3(cellnum,'VoxelIdxList');
    
    for hyb = 1:numChannels % hyb numbers
        for ch = 1:2 % fluorescent channel 1 or 2 for IF
            %mean_int = regionprops3(cellnum, Ipaint{hyb,ch}, 'MeanIntensity');
            
            if isempty(immunokey)
                target = '';
            else
                target = immunokey{(hyb-1)*2+ch}; % max 2 channels
            end
            
            for cellID = 1:size(s,1) % for each nucleus
                %fprintf(fileID,'%.0f,%.0f,%.0f,%s,%.0f,%.3f\n',...
                %        position, hyb, ch, target, cellID, mean_int.MeanIntensity(cellID));
                fprintf(fileID,'%.0f,%.0f,%.0f,%s,%.0f,%.3f,%.3f,%.3f,%.3f\n',...
                        position, hyb, ch, target, cellID, mean(Ipaint{hyb,ch}(s.VoxelIdxList{cellID})),...
                        quantile(Ipaint{hyb,ch}(s.VoxelIdxList{cellID}),0.25),...
                        median(Ipaint{hyb,ch}(s.VoxelIdxList{cellID})),...
                        quantile(Ipaint{hyb,ch}(s.VoxelIdxList{cellID}),0.75));
            end
        end
    end

    fclose(fileID);

end