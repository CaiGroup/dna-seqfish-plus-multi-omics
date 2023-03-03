function offsets = alignpointswrapper_v3lncrna(chArray, pointsPath, offsetsPath, chaTform, numRounds, folderArray, physicalTforms, usechabboffsets, usephysicaloffsets)
    % add options to apply chatforms or physical shift
    % tforms....automatically uses csv files to form tforms with is already
    % corrected for chabb corrections
    % align to initial fiducial marker instead of hyb1 RNA as initial fiducial
    % marker images were aligned to hyb1 DNA.

    %folderArray = 0:79; 
    %numRounds = 5; 
    offsets = cell(length(chArray), 1);
    pointsch = cell(length(chArray), 1);

    offsetsT = readtable(offsetsPath);
    for c = chArray
        %offsetsDir = fullfile(offsetsBaseDir, ['ch' num2str(c)], ['ch' num2str(c) '_offsets']);
        %offsetsPath = getfile(offsetsDir, ['offsets_ch' num2str(c)], 'match');
        %offsetsPath = getfile(pointsDir, ['offsets_ch' num2str(c)], 'match');
        offsets{c} = offsetsT(offsetsT.ch == c,:);
    end
    
    % set offsets equal to the first hyb, and align to ch1. as chromatic
    % aberration and physical shift are fixed and aligned to ch1.
    
    % x is the same as previous definition of col. so covert to it here.
    % -y is the same as previous definition of row. so covert to it here.
    
    % this should correct the offset unless all channels are failed.
    for h = folderArray+1
        if isnan(offsets{1,1}.x(h)) && ~isnan(offsets{2,1}.x(h))
            offsets{1,1}.x(h) = offsets{2,1}.x(h);
            offsets{1,1}.y(h) = offsets{2,1}.y(h);
            offsets{1,1}.z(h) = offsets{2,1}.z(h);
            
        elseif length(offsets)>2 && size(offsets{3},1)>0 && isnan(offsets{1,1}.x(h)) && ~isnan(offsets{3,1}.x(h))
            offsets{1,1}.x(h) = offsets{3,1}.x(h);
            offsets{1,1}.y(h) = offsets{3,1}.y(h);
            offsets{1,1}.z(h) = offsets{3,1}.z(h);
        end
        
        if isnan(offsets{2,1}.x(h)) && ~isnan(offsets{1,1}.x(h))
            offsets{2,1}.x(h) = offsets{1,1}.x(h);
            offsets{2,1}.y(h) = offsets{1,1}.y(h);
            offsets{2,1}.z(h) = offsets{1,1}.z(h);
        end
        
        if length(offsets)>2 && size(offsets{3},1)>0 && isnan(offsets{3,1}.x(h)) && ~isnan(offsets{2,1}.x(h))
            offsets{3,1}.x(h) = offsets{2,1}.x(h);
            offsets{3,1}.y(h) = offsets{2,1}.y(h);
            offsets{3,1}.z(h) = offsets{2,1}.z(h);
        end
        
    end      
    
    if length(offsets) == 2 || size(offsets{1},1) ~= size(offsets{3},1)% make third channel offset with second channel offset if offsets were calculated with two channels
        offsets{3,1} = offsets{2,1};
        offsets{3,1} = offsets{2,1};
        offsets{3,1} = offsets{2,1};
    end
    
    %ref_row = -offsets{1}.y(1);
    %ref_col = offsets{1}.x(1);
    %ref_z = offsets{1}.z(1);
    
    for c = chArray
        offsets{c}.row(:) = -offsets{c}.y(:);
        offsets{c}.col(:) = offsets{c}.x(:);
        offsets{c}.z(:) = offsets{c}.z(:);
    end

end

% combine this piece of code to the matlab function fro
% testNormalizeIntenisty as the full pipeline to get the points, aligne the
% images using the bead images, adjust the threshold ane decode the images