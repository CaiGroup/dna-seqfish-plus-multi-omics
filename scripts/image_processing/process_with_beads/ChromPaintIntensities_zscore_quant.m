function  [Fa2, Va2] = ChromPaintIntensities_zscore_quant(cellnum,Ipaint)
% function calculate Z-scored paint intensity per channel per cell
% optimized for 100k DNA seqFISH+ analysis
% Input should be 
% 1. cellnum, which is a 3D labeled image (same size as Ipaint image)
% 2. Ipaint (hyb61-96 chromosome paint images for 3 fluorescent channel)
% Date: 12/05/2020
% Author: Yodai Takei

numHybs = size(Ipaint,1);
numCh = size(Ipaint,2);

%IpaintZscore = cell(numHybs,numCh);
Fa2 = cell(numHybs,numCh);
Va2 = cell(numHybs,numCh);

for row = 1:numHybs
    for col = 1:numCh
        %indices_all = [];
        %Iz_all = [];
        % default image
        %Ib = zeros(size(Ipaint{row,col},1),size(Ipaint{row,col},2),size(Ipaint{row,col},3));
        for i = 1:max(max(max(cellnum))) % compute z-score label by label (nucleus by nucleus).
            indices  = find(cellnum==i); % grab a label(or nucleus) indices.
            I1 = Ipaint{row,col}(indices); % pull intensity in a label.
            Iz = zscore(double(I1)); % z-scored intensity values sorted by the indices in a label.
            %indices_all = [indices_all;indices]; % add indices
            %Iz_all = [Iz_all;Iz]; % add z-score
            
            Fa2{row,col}(i) = length(Iz(Iz>2))/length(Iz); % percentage of voxels z-score above 2 per marker per cell
            Va2{row,col}{i} = indices(Iz>2); % voxel indices z-score above 2 per marker per cell
        end
        %Ib(indices_all) = Iz_all; % map back the z-score (sorted by indices) to the image (Ib)
        %IpaintZscore{row,col} = Ib; % update IpaintZscore with the Z-scored image
    end
end

