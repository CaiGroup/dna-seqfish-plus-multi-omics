function L = roi2cellmask(vertex, imageSize)
% Makes a labeled mask from vertices of 2d ROIs
% 
% Author: Nico Pierson
% Date: October 11, 2018
% imageSize is usually [2048 2048 1] for a single .tif image

    % Set size cells
    numCells = length(vertex);
    L = zeros(imageSize);
    x = imageSize(1);
    y = imageSize(2);
    if length(imageSize) == 3
        z = imageSize(3);
    else
        z = 1;
    end

    for i = 1:numCells

        pm = poly2mask(vertex(i).x,vertex(i).y,x,y);
        pm3 = repmat(pm,[1 1 z]);

        % Add polymask to binaryMask
        L(pm3) = i; 
    end
    
    L = uint16(L);

end