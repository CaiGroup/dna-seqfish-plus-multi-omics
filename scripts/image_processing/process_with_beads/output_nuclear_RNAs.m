function output_nuclear_RNAs(experimentDir, nucsegDir, experimentLabel, position, ...
        alloweddiff, sqrtradius, channels)
    
    load(fullfile(nucsegDir,'segmentation', ['pos' num2str(position)], 'cellnum.mat') ,'cellnum');
    
    for channel = channels
        
        if (channel == 1) || (channel == 2)
            explabel = [num2str(alloweddiff) 'error-sqrt' num2str(sqrtradius) '-ch' num2str(channel)];
        else % channel 3
            explabel = ['ch' num2str(channel)];
        end
            
        savePath = fullfile(experimentDir, 'analysis', experimentLabel,  explabel);
        newPath = fullfile(savePath, 'NucCyto_RNAs');
        mkdir(newPath);
        listing = dir(fullfile(savePath,['*pos' num2str(position) '*.csv']));
        dotPath = fullfile(listing(1).folder, listing(1).name);
        points = readtable(dotPath);
        
        listSavePath = fullfile(newPath,[listing(1).name(1:end-4) '-NucCyto.csv']);
        fileID = fopen(listSavePath,'w');
        
        if (channel == 1) || (channel == 2) % barcoding     
            fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
            'fovID', 'cellID', 'geneID', 'regionID', 'x', 'y', 'z', 'seeds', 'intensity','NucCytoID'); % NucCytoID: nuclear RNA 0, cytoplasmic RNA 1
            TableName = ["cellID", "geneID", "regionID", "x", "y", "z", "seeds", "intensity"]; % in case header is not loaded properly.  
            points.Properties.VariableNames = TableName;
            for i = 1:size(points,1)
                nucID = cellnum(round(points.y(i)),round(points.x(i)),round(points.z(i)));
                if nucID == points.cellID(i) % filter out inconsist segmentation label dots (should be minimal or none)
                    fprintf(fileID,'%.0f,%.0f,%s,%.0f,%.3f,%.3f,%.3f,%.0f,%.3f,%.0f\n',...
                        position, points.cellID(i), convertCharsToStrings(points.geneID{i}), points.regionID(i), ...
                        points.x(i), points.y(i), points.z(i), points.seeds(i),points.intensity(i), 0);
                elseif nucID == 0 % dots outside nuclear segmentation
                    fprintf(fileID,'%.0f,%.0f,%s,%.0f,%.3f,%.3f,%.3f,%.0f,%.3f,%.0f\n',...
                        position, points.cellID(i), convertCharsToStrings(points.geneID{i}), points.regionID(i), ...
                        points.x(i), points.y(i), points.z(i), points.seeds(i),points.intensity(i), 1);
                end
            end
            
        else % sequential
            fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
                    'fovID', 'cellID', 'hybID', 'geneID', 'x', 'y', 'z', 'intensity','NucCytoID');   
            for i = 1:size(points,1)
                nucID = cellnum(round(points.y(i)),round(points.x(i)),round(points.z(i)));
                if nucID == points.cellID(i)
                    fprintf(fileID,'%.0f,%.0f,%.0f,%s,%.3f,%.3f,%.3f,%.3f,%.0f\n',...
                        points.fov(i), points.cellID(i), points.hybID(i), convertCharsToStrings(points.geneID{i}), ...
                        points.x(i), points.y(i), points.z(i), points.intensity(i), 0);
                elseif nucID == 0
                    fprintf(fileID,'%.0f,%.0f,%.0f,%s,%.3f,%.3f,%.3f,%.3f,%.0f\n',...
                        points.fov(i), points.cellID(i), points.hybID(i), convertCharsToStrings(points.geneID{i}), ...
                        points.x(i), points.y(i), points.z(i), points.intensity(i), 1);
                end
            end
            
        end
        
        fclose(fileID);
        
    end

end