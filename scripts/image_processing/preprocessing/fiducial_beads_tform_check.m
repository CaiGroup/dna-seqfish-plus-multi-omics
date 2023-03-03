function I = fiducial_beads_tform_check(experimentDir, position, finalFixedPath)
% 
% Date: 09/08/2021
% Author: Yodai Takei (modified original script from Nico Pierson)
% ytakei@caltech.edu

    tic
    %% Check if Fiji and bfmatlab is in the path
    fijiDirectory = checkfijipath();
    checkbfmatlabpath();
    
    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    endingDateString = datestr(dateStart, formatDate);
    
    
    %% Initialize Variables
    folderSaveName = 'processedimages_rotalign';


    % Make save directory: default 'organizehybs' folder
    saveDir = fullfile(experimentDir, folderSaveName, ['pos' num2str(position)]);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    
    
    %% Get the tform for RNA to DNA experiment
    imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
    imFixedPath = fullfile(finalFixedPath, imageName);
    % initial RNA fiducial marker image without DAPI (3 channels)
    [I, sizeCDNA, sizeZDNA, ~, ~] = grabimseries(imFixedPath, position);
    
    toc
end