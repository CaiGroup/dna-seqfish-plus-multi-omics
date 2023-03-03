function L = getlabel(filepath, varargin)
% returns labeled image from roi(.zip) or h5(.h5) file
%
% filepath is the path of the roi.zip or the .h5 labeled image
% z is the number of z-slices
% date: 2/21/2019

    %% Set up optional Parameters
    argsLimit = 1;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('seqfish:getlabel:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    % set defaults for optional inputs
    optargs = {1};
    optargs(1:numvarargs) = varargin;
    [z] = optargs{:};
    
    %% get the label
    imagesize = [2048 2048 z];
    ext = filepath(end-2:end);
    if strcmp(ext,'zip')
        vertex = selfsegzip(filepath);
        L = roi2cellmask(vertex, imagesize);
    elseif strcmp(ext,'.h5')
        L = geth5mask(filepath);
        %{
    elseif strcmp(ext,'roi')
        % if multiple roi, then break otherwise only one roi
        warning 'wrong file if more than one .roi\n';
        vertex = selfseg(filepath);
        L = roi2cellmask(vertex, imagesize);
        %}
    elseif strcmp(ext,'mat')
        t = load(filepath);
        L = t.cellnum;
        % adjust the label: rep3-2, but not rep2-2
        % label 0 and 1 is the background, set 2 to label 1
        L = L - 1;
        L(L<0) = 0;
    end

end