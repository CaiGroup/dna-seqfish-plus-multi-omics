function Ifoci = foci_preprocess(I)

% modified version for 100k IF.

%addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts');

%saveDir = fullfile(experimentDir, 'immuno-data', 'focinum-r9');
%mkdir(saveDir);
%saveDir2 = fullfile(experimentDir, 'immuno-data', 'focinum-r9',['pos' num2str(position)]);
%mkdir(saveDir2);
%saveImName = ['focinum-Pos' num2str(position) '-hyb' num2str(min(cycleArray)) '-' num2str(max(cycleArray)) '.mat'];
%saveImPath = fullfile(saveDir, saveImName);

numCh = size(I,2);
numHybs = size(I,1);
Ifoci = I;
%key_num = 1;
for row = 1:numHybs
    for col = 1:numCh
        Miji;
        MIJ.createImage('Ib',I{row,col},true);
        MIJ.run('Subtract Background...', 'rolling=9 stack');
        MIJ.run('Gaussian Blur...', 'sigma=1 stack');
        MIJ.run('Auto Threshold', 'method=Yen white stack');
        %MIJ.run('Dilate', 'stack');
        %MIJ.run('Erode', 'stack');
        Ib = uint16(MIJ.getCurrentImage) >0;
        MIJ.run('16-bit');
        %MIJ.createImage('Ib2',I{row,col},true);
        %MIJ.run('Merge Channels...', 'c1=Ib2 c2=Ib create');
        %MIJ.run('Save', ['save=[' saveDir2 '\pos' num2str(position) '-' num2str(key_num) '-' immunokey{key_num} '-foci.tif' ']']);
        MIJ.run('Close All')
        MIJ.exit
        %L = bwareaopen(Ib,1); % lower threshold for the voxel size
        L = bwareaopen(Ib,20); % 20 for the r = 9 condition
        L = L-bwareaopen(L,10000000); % upper threshold for the voxel size
        focinum = bwlabeln(L);
        Ifoci{row,col} = focinum;
        %[y,x,z] = ind2sub(size(focinum),find(focinum ~= 0));
        %linInd = sub2ind(size(focinum),y,x,z);
        %fociL = focinum(linInd);
        %key_num = key_num+1;
    end
end








