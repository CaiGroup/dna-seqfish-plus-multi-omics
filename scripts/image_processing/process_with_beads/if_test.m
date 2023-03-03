function if_test(experimentDir, experimentLabel, position, numChannels, IpaintZscore, save_folder, save_name)

listSavePath = fullfile(experimentDir,'analysis',experimentLabel,'mean-IF-intensity',['Zscore-test3-pos' num2str(position) '.csv']);
fileID = fopen(listSavePath,'w');
if_num = 1:numChannels*2; % total target
header_name = ["fov" "hyb" "channel" "cellID" "x" "y" "z" "dot_int" "Round 1"...
       "Round 2" "Round 3" "Round 4" "name" "chrom"];
for t = 1:length(if_num)
    header_name = [header_name "IF_"+string(t)];
end

% for header
fmt=[repmat('%s,',1,length(header_name)) '\n'];
% for actual datasets
fmt2=['%.0f,%.0f,%.0f,%.0f,%.3f,%.3f,%.3f,%.0f,%.0f,%.0f,%.0f,%.0f,%s,%s,' repmat('%3f,',1,length(if_num)) '\n'];

fprintf(fileID, fmt, header_name);

decoded_Path = fullfile(save_folder, [save_name num2str(position) '.csv']);
decoded_points = readtable(decoded_Path);

index = decoded_points.fov~=position;
decoded_points(index,:) = [];

for dot = 1:size(decoded_points,1)
   
   v1 =   decoded_points.hyb(dot);  
   v2 =   decoded_points.channel(dot);  
   v3 =   decoded_points.cellID(dot);  
   v4 =   decoded_points.dot_int(dot);  
   v5 =   decoded_points.Round1(dot);  
   v6 =   decoded_points.Round2(dot);  
   v7 =   decoded_points.Round3(dot);  
   v8 =   decoded_points.Round4(dot);  
   v9 =   string(decoded_points.name(dot));  
   v10 =  string(decoded_points.chrom(dot));  
            
   x =   decoded_points.x(dot);  
   y =   decoded_points.y(dot);  
   z =   decoded_points.z(dot);  
   
   x_r = round(x);
   y_r = round(y);
   z_r = round(z);
   
   header_name = ["fov" "hyb" "channel" "cellID" "x" "y" "z" "dot_int" "Round 1"...
       "Round 2" "Round 3" "Round 4" "name" "chrom"];
   
   out_list = [position, v1, v2, v3, x, y, z, v4, v5, v6, v7, v8, v9, v10]; %make as string list?
   out_list = string(out_list);
   
   for hyb = 1:numChannels
       for ch = 1:2
           if_int = IpaintZscore{hyb,ch}(y_r,x_r,z_r);
           out_list = [out_list if_int];
       end
   end
   
   fprintf(fileID,fmt2,out_list);

end

fclose(fileID);
