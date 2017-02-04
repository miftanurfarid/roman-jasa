function generate_signals(X, Az, output_dir)

%Simulate acoustic mixture for the two-source scenario with locations at az and given input files.
%Output: Binaural signals for the individual sources and for the mixture.


nrsource=length(Az);
if size(X,1)~=nrsource
    'eror: X(nrsource,max_signal)'
end

name_l=strcat(output_dir,'\Left');
name_r=strcat(output_dir,'\Right');
name_l1=strcat(output_dir,'\Left1');
name_r1=strcat(output_dir,'\Right1');
name_l2=strcat(output_dir,'\Left2');
name_r2=strcat(output_dir,'\Right2');


addpath HRTF\matlab_scripts\;

% signal coming from az1

az=Az(1);
[hrtf]=readhrtf(0,abs(az),'H');

if az>=0
hl=[hrtf(1,:)];
hr=[hrtf(2,:)];
else
hl=[hrtf(2,:)];
hr=[hrtf(1,:)];
end

x=X(1,:);
left1=conv(x,hl);
right1=conv(x,hr);

    
left2=zeros(size(left1));
right2=zeros(size(right1));
   
% interference

for n=2:nrsource

az=Az(n);
[hrtf]=readhrtf(0,abs(az),'H');

if az>=0
hl=[hrtf(1,:)];
hr=[hrtf(2,:)];
else
hl=[hrtf(2,:)];
hr=[hrtf(1,:)];
end

x=X(n,:);
left2=left2+conv(x,hl);
right2=right2+conv(x,hr);

end

% the sum of the 2 signals 


left=left1+left2;

right=right1+ right2;


f=fopen( name_l,'w');
fprintf(f,'%f\n ',left);
fclose(f);

f=fopen(name_r,'w');
fprintf(f,'%f\n',right);
fclose(f);

if nrsource>1

f=fopen( name_l1,'w');
fprintf(f,'%f\n ',left1);
fclose(f);

f=fopen(name_r1,'w');
fprintf(f,'%f\n',right1);
fclose(f);

f=fopen( name_l2,'w');
fprintf(f,'%f\n ',left2);
fclose(f);

f=fopen(name_r2,'w');
fprintf(f,'%f\n',right2);
fclose(f);

end





