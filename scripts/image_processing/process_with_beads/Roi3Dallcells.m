function BW = Roi3Dallcells(segPath, position)
    listing = dir([segPath '\Pos' num2str(position) '\*cellnum*.mat']);
    load([listing(length(listing)).folder '\' listing(length(listing)).name],'cellnum');
    % assuming last file is the latest with the time stamps.
    BW = zeros(size(cellnum));
    indices = find(cellnum>0);
    BW(indices) = 1; 
end