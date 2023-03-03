function [offsets] = alignpointswrapper_v4IF(chArray, pointsPath, offsetsPath, chaTform, numRounds, folderArray, physicalTforms, usechabboffsets, usephysicaloffsets, offsetsPath_ref, hyb1_ref)
    % add options to apply chatforms or physical shift --> shouldn't be
    % used in the latest version.
  

    %folderArray = 0:79; 
    %numRounds = 5; 
    offsets = cell(length(chArray), 1);
    offsets_ref = cell(length(chArray), 1);
    pointsch = cell(length(chArray), 1);

    offsetsT = readtable(offsetsPath);
    for c = chArray
        %offsetsDir = fullfile(offsetsBaseDir, ['ch' num2str(c)], ['ch' num2str(c) '_offsets']);
        %offsetsPath = getfile(offsetsDir, ['offsets_ch' num2str(c)], 'match');
        %offsetsPath = getfile(pointsDir, ['offsets_ch' num2str(c)], 'match');
        offsets{c} = offsetsT(offsetsT.ch == c,:);
    end
    
    offsetsT_ref = readtable(offsetsPath_ref);
    for c = chArray
        offsets_ref{c} = offsetsT_ref(offsetsT_ref.ch == c,:);
    end
    
    % For IF, the alignment was done only with ch3 so copy it to ch1,2
    offsets{1,1} = offsets{3,1};
    offsets{2,1} = offsets{3,1};
    
    % set offsets equal to the first hyb, and align to ch1. as chromatic
    % aberration and physical shift are fixed and aligned to ch1.
    
    % x is the same as previous definition of col. so covert to it here.
    % -y is the same as previous definition of row. so covert to it here.
        
    if hyb1_ref
        % align to hyb1 of DNA seqFISH    
        for c = chArray
            offsets{c}.row(:) = -offsets{c}.y(:) + offsets_ref{3}.y(1) ; % align to hyb1 of ch3 as chromatic shift is already corrected in the latest version.
            offsets{c}.col(:) = offsets{c}.x(:) - offsets_ref{3}.x(1);
            offsets{c}.z(:) = offsets{c}.z(:) - offsets_ref{3}.z(1);
        end
    else
        % align to initial_fiducial_markers of DNA seqFISH or RNA seqFISH   
        for c = chArray
            offsets{c}.row(:) = -offsets{c}.y(:);
            offsets{c}.col(:) = offsets{c}.x(:);
            offsets{c}.z(:) = offsets{c}.z(:);
        end
    end

end

% combine this piece of code to the matlab function fro
% testNormalizeIntenisty as the full pipeline to get the points, aligne the
% images using the bead images, adjust the threshold ane decode the images