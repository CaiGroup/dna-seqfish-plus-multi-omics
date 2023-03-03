function finalpoints_paints = ChromPaintIntensities_perpoints_v4(experimentDir, experimentName, experimentLabel, finalpoints, I, position, segment)
% This works for ch3 analysis for DNA seqFISH+ experiments.
% Changed format to pass variables from seqfish pipeline.
% points: (1:60,1) hyb1-60 points (already aligned and corrected for shifts)
% I: (1:20,1) hyb61-80 of ch3 for chromosome paint. (already aligned and corrected for shifts)
% integrated with roi = 'whole' for the 3D ROI analysis.
%
% Author: Yodai Takei
% Date updated: 04/08/20
% Email: ytakei@caltech.edu

switch segment
    case 'roi'
        vertex = selfsegzip([experimentDir, '\segmentation_2Dtest\pos' num2str(position) '\RoiSet.zip']);
        numCells = length(vertex);
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

        if numCells > 1 % only for the 2D ROIs
            for i =1:numCells
            polygonIn = inpolygon(finalpoints_paints{hyb,channel}.channels(:,1),finalpoints_paints{hyb,channel}.channels(:,2),vertex(i).x,vertex(i).y);
                for dot = 1:length(polygonIn)
                    if polygonIn(dot) ~= 0
                        finalpoints_paints{hyb,channel}.cellID(dot,1) = i*polygonIn(dot);
                    end
                end
            end
        end
    end

    saveDirPath = fullfile(experimentDir, 'analysis', experimentLabel, 'paint decoding');
    if exist(saveDirPath, 'dir') ~= 7
        mkdir(saveDirPath);
    end

    savePath = fullfile(saveDirPath, ['points-paints-pos' num2str(position) '-' experimentName '.mat']);
    save(savePath, 'finalpoints_paints', 'position');
end


