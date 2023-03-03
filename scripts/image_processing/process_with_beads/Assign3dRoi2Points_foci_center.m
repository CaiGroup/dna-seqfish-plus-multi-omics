function [cell_id, Stats] = Assign3dRoi2Points_foci_center(cellnum,finalPosList)

% to add an additional column of cellnum to input from .csv file.
% finalPosList should be table format.

gene_num = size(finalPosList,1);
%name_gene = [];
Stats = [];
stats = regionprops3(cellnum);
Stats(:,4) = stats.Volume;
Stats(:,1) = stats.Centroid(:,1);
Stats(:,2) = stats.Centroid(:,2);
Stats(:,3) = stats.Centroid(:,3);
% 1 is already filtered out with getlabel function. so commented this.
%Stats(1,:) = [];% 1 is outside cell ROI, Column 1,2,3,4 represent x,y,z centroid and volume (voxel)


points_gene_original = finalPosList.Centroid;
points_gene = round(points_gene_original);
indices = find(points_gene(:,1)>size(cellnum,1)|points_gene(:,1)<1|points_gene(:,2)>size(cellnum,2)|points_gene(:,2)<1|points_gene(:,3)>size(cellnum,3)|points_gene(:,3)<1);
points_gene(indices,:) = [];
linInd = sub2ind(size(cellnum),points_gene(:,2),points_gene(:,1),points_gene(:,3));
cells_gene = cellnum(linInd); %returned the cell ID, starting from 2
filter_indices = 1;
cell_id = zeros(gene_num,1);
for i = 1:gene_num
    if isempty(find(indices==i, 1))
        cell_id(i) = cells_gene(filter_indices); %cell ID starts from 1. 0 is outside. already finalized indices. so no need to subtract 1.
        filter_indices = filter_indices + 1;
    end
end

%savePathCh = fullfile('F:\Yodai\DNA+\2019-09-09-brain-rep2-2-DNAFISH', 'brain-rep2-2-DNAFISHdecoded-ch2-Pos1-20191026.mat');
%save(savePathCh, 'cell_id', 'Stats', 'cellnum_filename');

