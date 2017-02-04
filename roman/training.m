function training(namefile_ITD, namefile_IID, namefile_Ratio, Az, namefile_configuration)

%Generate MAP decision rules for all channels based on the training data
%Input:
%  TRAINFILE_ITD, TRAINFILE_IID, TRAINFILE_R = estimates of ITD, IID and the a priori relative strength for all channels based on the training data.
%Output:
%  D = surface decisions
%  Region = feature domain information

'training'

Npoints=100;
display=0;
max_channel=128;

D=zeros(Npoints,Npoints,max_channel);
Region=zeros(5,max_channel);


P=load(namefile_ITD);
u=rand(size(P));x=norminv(u,0,0.3);P=P+x;
P=reshape(P,max_channel,length(P)/max_channel);
L=load(namefile_IID);L=reshape(L,max_channel,length(L)/max_channel);
R=load(namefile_Ratio);R=reshape(R,max_channel,length(R)/max_channel);


for chan=1:max_channel
    
chan

n=size(R,2);
I1=find(R(chan,:)>0.5);
I2=find(R(chan,:)<=0.5);
len1=length(I1);len2=length(I2);

if len1>2000
      len1=2000;
end
if len2>2000
      len2=2000;
end
             
l1=round(rand(1,len1)*(n-1)+1);
l2=round(rand(1,len2)*(n-1)+1);     
  
I=[l1,l2];
length(l1)
length(l2)

features=[P(chan,I)-45;L(chan,I)];
targets=(R(chan,I)>0.5);

size(find(targets==1))
size(find(targets==0))

region=[min(features(1,:)),max(features(1,:)),min(features(2,:)),max(features(2,:)),Npoints];
d=decisionRegion(features,targets,region,display);
D(:,:,chan)=d;
Region(:,chan)=region';

end


save (namefile_configuration, 'D', 'Region','Az')