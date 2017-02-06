function estimateMask(input_directory, D, Region)
% Compute estimated binary mask from from itd/iid features using the learned decisionrules

N = Region(5,1);

if isunix == 1
    dir_ch = '/';
end

if ispc == 1
    dir_ch = '\';
end

inputfile_itd = strcat(input_directory,dir_ch,'itd');
inputfile_iid = strcat(input_directory,dir_ch,'iid');
outputfile_mask = strcat(input_directory,dir_ch,'mask');


iid = load(inputfile_iid);
iid = reshape(iid,128,length(iid)/128);

itd = load(inputfile_itd);
itd = reshape(itd,128,length(itd)/128);

fprintf('size iid = %s x %d\n',size(iid))

l = size(iid,2);
max_channel = size(iid,1);

itd = itd-45;

mask = zeros(size(itd));
        
for i = 1:l
    for j = 1:max_channel
        L = round(1/(Region(4,j)-Region(3,j))*...
            ((N-1)*iid(j,i)+Region(4,j)-N*Region(3,j)));
        
        P = round(1/(Region(2,j)-Region(1,j))*...
            ((N-1)*itd(j,i)+Region(2,j)-N*Region(1,j)));
        
        P = min(max(1,P),N);
        
        L = min(max(1,L),N);
      
        mask(j,i) = D(L,P,j);
    end
end
 
f = fopen(outputfile_mask,'w');
fprintf(f,'%d\n',l);

for i = 1:l
    for j = 1:max_channel
        fprintf(f,'%f ',mask(j,i));
    end
    fprintf(f,'\n',1);
end

fclose(f);

end

% inputfile_ratio=strcat(input_directory,'/ratio');
% ratio=load(inputfile_ratio);
% ratio=reshape(ratio,128,length(ratio)/128);
% 
% figure;imagesc(abs(mask)>0)