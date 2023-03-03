function getifv2_Va2v2(position, cellnum, Va2, numChannels, experimentDir, experimentLabel, experimentDir_ref2)
% function to output mean and quantile IF intensity per marker per cell as a csv file
% cellnum: labeled image
% Ipaint: aligned IF images
% numChannels: total hyb number
% Immunokey file should be in "immunokey" folder and sorted by hyb and channel
% and the first column should be antibody names.
% Date: 12/17/2021
% Author: Yodai Takei
    
    saveGlobalDir = fullfile(experimentDir,'analysis',experimentLabel,'IF-Va2-per-cell');
    if exist(saveGlobalDir, 'dir') ~= 7
        mkdir(saveGlobalDir);
    end
    
    % load sequential codebook (IF + lncRNA)
    immunoPath = getfile(fullfile(experimentDir_ref2, 'immunokey'), '.', 'match');
    if exist(immunoPath, 'file') ~= 2
        immunokey = [];
    else
        immunokey = readsequential(immunoPath, 'header');
    end
    
    listSavePath = fullfile(saveGlobalDir, ['output-IF-Va2-per-cell-pos' num2str(position) '.csv']);
    fileID = fopen(listSavePath,'w');
    fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
                'fov', 'hyb_1','hyb_2','channel_1','channel_2','target_1','target_2','cellID', 'percentage');
            
    s = regionprops3(cellnum,'VoxelIdxList');
    
    % compare markers between channels 1 and 2
    for hyb = 1:numChannels % hyb numbers
        for hyb2 = 1:numChannels % hyb numbers

            if isempty(immunokey)
                target1 = '';
                target2 = '';
            else
                target1 = immunokey{(hyb-1)*2+1}; % max 2 channels
                target2 = immunokey{(hyb2-1)*2+2}; % max 2 channels
            end

            for cellID = 1:size(s,1) % for each nucleus

                v1 = Va2{hyb,1}{cellID};
                v2 = Va2{hyb2,2}{cellID};
                
                Lin = length(intersect(v1,v2));
                Lall = length(unique([v1;v2]));
                p = Lin/Lall*100;

                fprintf(fileID,'%.0f,%.0f,%.0f,%.0f,%.0f,%s,%s,%.0f,%.3f\n',...
                        position, hyb, hyb2, 1, 2, target1, target2, cellID, p);
            end
        end
    end
    
    % compare markers within channels
    for hyb = 1:numChannels % hyb numbers
        for hyb2 = hyb:numChannels % hyb numbers
            for ch = 1:2
                if isempty(immunokey)
                    target1 = '';
                    target2 = '';
                else
                    target1 = immunokey{(hyb-1)*2+ch}; % max 2 channels
                    target2 = immunokey{(hyb2-1)*2+ch}; % max 2 channels
                end

                for cellID = 1:size(s,1) % for each nucleus

                    v1 = Va2{hyb,ch}{cellID};
                    v2 = Va2{hyb2,ch}{cellID};

                    Lin = length(intersect(v1,v2));
                    Lall = length(unique([v1;v2]));
                    p = Lin/Lall*100;

                    fprintf(fileID,'%.0f,%.0f,%.0f,%.0f,%.0f,%s,%s,%.0f,%.3f\n',...
                            position, hyb, hyb2, ch, ch, target1, target2, cellID, p);
                end
            
            end
        end
    end

    fclose(fileID);

end