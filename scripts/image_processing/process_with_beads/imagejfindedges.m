function processImage = imagejfindedges(image, varargin)
% use the imageJ Find Edges to get the edges
%
% requirements:
% add the bfmatlab package bioformats_package.jar to the matlab
% javaclasspath
% bfmatlab = 'C:\Users\Long Cai - 2\Box\code\io\bfmatlab\bioformats_package.jar';
% javaaddpath(bfmatlab);
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 5/19/2020
% modified the output format by YT on 01/19/22

    %% Set up optional Parameters
    
    numvarargs = length(varargin);
    if numvarargs > 2
        error('myfuns:getallbeads:TooManyInputs', ...
            'requires at most 2 optional inputs');
    end

    % set defaults for optional inputs
    optargs = {'tempProcessed', []}; % default of using 7 x 7 pixel grid for gaussian function
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Place optional args in memorable variable names
    [uniqueString, fijiDirectory] = optargs{:};

    %% Declare Variables
    % instead of using fijiDirectory, use temporary directory
    if isempty(fijiDirectory)
        fijiDirectory = pwd;
    end
    %tempDirectory = tempdir; % temp direcotyr
    
    % Start Fiji in the background
    Miji(false); 
    
    %% Create image in ImageJ
    namesh = ['C' num2str(1) '-'  num2str(1) '.tif'];
    MIJ.createImage(namesh, image, true);
    % add .tif to uniqueString
    uniqueString = [uniqueString '.tif'];
    saveBackSubImPath = fullfile(fijiDirectory, uniqueString);

    %% Subtract Background using Rolling Ball Algorithm in ImageJ
    try
        MIJ.run("Find Edges", "stack");
        
        MIJ.run('Save', ['save=[' saveBackSubImPath ']']); % use a unique string to temporarily save image
        MIJ.run('Close All')
    catch
        error('MIJ exited incorrectly: most likely caused by out of memory in the java heap\n');
    end
    
    tempPosition = 0;
    [processImage, ~, ~, ~, ~] = grabimseries(fullfile(fijiDirectory, uniqueString), tempPosition);
    % try to delete temp file
    warning('off','all');
    delete(fullfile(fijiDirectory, uniqueString));
    warning('on','all');
    processImage = processImage{1,1};

end


