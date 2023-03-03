function I = applydapitform(image, tform)    
% applys tform to image in cell array
%
% Return image with applied transformation
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 2/26/2019
% Modified:

    % Variables
    numChannels = length(image);
    I = cell(1, numChannels);
    numDim = length(tform.T);
    if numDim ~= 3 && numDim ~= 4
        error 'tform does not have the right number of dimensions';
    end

    if numDim == 4
        for i = 1:numChannels
                I{i} = imwarp(image{i}, tform, 'OutputView', imref3d(size(image{i})));
        end
    elseif numDim == 3
        numZ = size(image{1},3);
        for i = 1:numChannels
            for z = 1:numZ
                I{i}(:,:,z) = imwarp(image{i}(:,:,z), tform, 'OutputView', imref2d(size(image{i}(:,:,z))));
            end
        end
    end
    
    
end