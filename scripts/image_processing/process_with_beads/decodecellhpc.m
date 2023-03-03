function [] = decodecellhpc(points, experimentName, position, channel, sqrtradius, ...
    alloweddiff, numRounds, numChannels, cell, saveDir, segmentPath, segment, barcodekey, iter, offpercent, minseeds)
% Wrapper function to load points and call function decodeallcells.m to
% decode each cell
%
% Structure:
% 1. use perl script to create loops - calls sh function
% 2. sh function creates a job by calling the matlab function to decode the
% data for each cell in each position
% 3. Need a function to combine everything togeter and calculate the false
% positive rate
%
% Need to call output function separately once the function is done
% processing everything.

        
        %% decode all the cells in the position for each cell
        % segment cells - for rois
        [pointsxcell, numCells, numpointspercell] = segmentpoints2cells(segmentPath, points, segment);

        % decode the points in the cell
        [dotlocations, seeds] = decodeimageshpc(pointsxcell{cell}, ...
            barcodekey, numRounds, numChannels, alloweddiff, sqrtradius);
        
        
        %% loop over min number of seeds and save as a separate file
        
        for minseed = minseeds
            %% Filter the Based on the number of Seeds 
            [finalPosList, PosList, dotlocations, numdotlocations, numpointconsensus, numfinalpoints] = filterseedsv2(seeds, dotlocations, minseed);


            % save the data
            explabel = [num2str(alloweddiff) 'error-sqrt' num2str(sqrtradius) '-iter' num2str(iter) '-ch' num2str(channel) '-.' num2str(offpercent) 'offpercent'];
            savePath = fullfile(saveDir, ['decodeData-minSeeds' num2str(minseed) '-Pos' num2str(position) '-Cell' num2str(cell) '-' explabel '-' experimentName '.mat']);
            save(savePath, 'finalPosList', 'posList', 'dotlocations', 'numpointconsensus', ...
                'numdotlocations', 'numfinalpoints', 'seeds', 'numpointspercell', ...
                'numCells', 'pointsxcell', 'points', 'position', 'cell', 'barcodekey', 'minseed');
        end
end




