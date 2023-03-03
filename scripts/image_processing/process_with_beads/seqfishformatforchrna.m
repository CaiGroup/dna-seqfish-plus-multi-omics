function points = seqfishformatforchrna(rawpoints, channel, numRounds, numChannels, folderArray)
% formatting the pointsch for mRNA/intron RNA seqFISH input
% Author: Yodai Takei
% Date updated: 09/13/21
% Email: ytakei@caltech.edu

points = cell(numRounds, 1);
for i = 1:numRounds
    points{i} = struct('channels', cell(numChannels,1), 'intensity', cell(numChannels, 1), 'scaledIntensity', cell(numChannels,1));
end

for folder = folderArray
    idx = folder + 1;
    r = ceil(idx / numChannels);
    ch = mod(idx, numChannels);
    if ch == 0
        ch = numChannels;
    end
    points{r}(ch).channels = rawpoints{channel,1}{1,1}(idx).channels;
    points{r}(ch).intensity = rawpoints{channel,1}{1,1}(idx).intensity/median(rawpoints{channel,1}{1,1}(idx).intensity);
    points{r}(ch).scaledIntensity = rawpoints{channel,1}{1,1}(idx).scaledIntensity;
end