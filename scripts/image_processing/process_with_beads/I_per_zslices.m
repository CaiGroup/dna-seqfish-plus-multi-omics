function Iout = I_per_zslices(I)

% format I (images) to new cell Iout z by z
% written by Yodai Takei 
% Last updated on 12/17/21.

numHybs = size(I,1);
numChannels = size(I,2);
numZ = size(I{1,1},3);

Iout = cell(numZ,1);

for z = 1:numZ
    Iout{z} = cell(numHybs,numChannels);
    for hyb = 1:numHybs
        for ch = 1:numChannels
            Iout{z}{hyb,ch} = I{hyb,ch}(:,:,z);
        end
    end
end