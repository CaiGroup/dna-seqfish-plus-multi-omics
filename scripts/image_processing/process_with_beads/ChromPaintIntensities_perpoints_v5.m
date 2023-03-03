function finalpoints_paints = ChromPaintIntensities_perpoints_v5(experimentDir, experimentName, experimentLabel, finalpoints, I, position, segment)
% This is for the 100k DNA seqFISH+ analysis to assign chromosome paint
% intensity profiles (hyb61-96) to each DNA seqFISH+ spots indentified in
% hyb1-60, and assign cellID with segmentation labels.
% Changed format to pass variables from seqfish pipeline.
% points: (1:60,1:3) hyb1-60 points (already aligned and corrected for shifts)
% I: (1:36,1:3) hyb61-96 for chromosome paint. (already aligned and corrected for shifts)
% updated segmentation option (segment) as 'roi':imageJ 2D ROI, 'cellpose':labeled 3D tif image, 'whole': without segmentation input.
% corresponding segmentation files should be in the folder defined below.
%
% Author: Yodai Takei
% Date updated: 06/23/21
% Email: ytakei@caltech.edu

switch segment
    case 'roi'
        vertex = selfsegzip([experimentDir, '\segmentation\pos' num2str(position) '\RoiSet.zip']);
        numCells = length(vertex);
     
    case 'cellpose'
        imBeadsPath = fullfile(experimentDir, 'segmentation', ['MMStack_Pos' num2str(position)], 'Segmentation', 'labeled_img_post.tif');
        [segImsTemp, ~, ~, ~, ~] = grabimseries(imBeadsPath, position);
        cellnum = double(segImsTemp{1,1});
        numCells = max(max(max(cellnum)));
        
    case 'cellnum'
        segPath = fullfile(experimentDir, 'segmentation', ['pos' num2str(position)], 'cellnum.mat');
        load(segPath,'cellnum');
        numCells = max(max(max(cellnum)));
        
    case 'whole'
        numCells = 1;
end


finalpoints_paints = finalpoints;
x_max = size(I{1},2);
y_max = size(I{1},1);
z_max = size(I{1},3);
indices = [];

for channel = 1:size(finalpoints,2)
    
    points = finalpoints(:,channel);

    for hyb = 1:length(points)
        idx = 1;
        for dot = 1:length(points{hyb,1}.channels)
            location = round(points{hyb,1}.channels(dot,:));
            if (points{hyb,1}.intensity(dot)~=0)&&(1<=location(2))&&(location(2)<=x_max)&&(1<=location(1))&&(location(1)<=y_max)&&(1<=location(3))&&(location(3)<=z_max) % to remove dots outside the images.
                for i = 1:size(I,1) % chromosome paints
                    finalpoints_paints{hyb,channel}.intensitypaint(dot,i) = I{i,channel}(location(2),location(1),location(3));
                end
            else
                indices(idx) = dot;
                idx = idx+1;
                for i = 1:size(I,1)
                    finalpoints_paints{hyb,channel}.intensitypaint(dot,i) = 0; % in case for last rows
                end
            end
        end


        % filter out problematic dots above.

        finalpoints_paints{hyb,channel}.channels(indices,:) = [];
        finalpoints_paints{hyb,channel}.intensity(indices,:) = [];
        finalpoints_paints{hyb,channel}.intensitypaint(indices,:) = [];

        indices = [];

        finalpoints_paints{hyb,channel}.cellID(:,1) = zeros(length(finalpoints_paints{hyb,channel}.channels),1);
        switch segment 
            case 'roi'
                if numCells > 1 
                    for i =1:numCells
                    polygonIn = inpolygon(finalpoints_paints{hyb,channel}.channels(:,1),finalpoints_paints{hyb,channel}.channels(:,2),vertex(i).x,vertex(i).y);
                        for dot = 1:length(polygonIn)
                            if polygonIn(dot) ~= 0
                                finalpoints_paints{hyb,channel}.cellID(dot,1) = i*polygonIn(dot);
                            end
                        end
                    end
                end
                
            case {'cellpose','cellnum'} % find cellpose label on a rounded voxel.
                if numCells > 1 
                    for dot = 1:length(finalpoints_paints{hyb,channel}.channels)
                        finalpoints_paints{hyb,channel}.cellID(dot,1) = cellnum(round(finalpoints_paints{hyb,channel}.channels(dot,2)),round(finalpoints_paints{hyb,channel}.channels(dot,1)),round(finalpoints_paints{hyb,channel}.channels(dot,3)));
                    end
                end
                
            case 'whole'
            % all numbers will be 0.
        end
    end
end

saveDirPath = fullfile(experimentDir, 'analysis', experimentLabel, 'paint decoding');
if exist(saveDirPath, 'dir') ~= 7
    mkdir(saveDirPath);
end

savePath = fullfile(saveDirPath, ['points-paints-pos' num2str(position) '-' experimentName '.mat']);
save(savePath, 'finalpoints_paints', 'position','-v7.3');

