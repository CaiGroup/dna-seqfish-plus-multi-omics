function finalpoints = seqfishformatforch123(pointsch, numpointshyb)
% formatting the pointsch to finalpoints
% Author: Yodai Takei
% Date updated: 12/07/20
% Email: ytakei@caltech.edu

numHybs = numpointshyb; % number of loci imaging rounds (= 60)
numChannels = length(pointsch); % number of fluorescent channel (= 3)
finalpoints = cell(numHybs,numChannels);

for ro = 1:numHybs
    for bch = 1:numChannels
        finalpoints{ro,bch}.channels = pointsch{bch,1}{1,1}(ro).channels;
        finalpoints{ro,bch}.intensity = pointsch{bch,1}{1,1}(ro).intensity;
    end
end