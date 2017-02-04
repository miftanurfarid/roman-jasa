
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%TESTING EXAMPLE - TWO-SOURCE CONFIGURATION%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_testing(namefile_configuration, inputfile_target, inputfile_noise1,inputfile_noise2)

 F1=441; %sampling frequency of the HRTF database is 44.1 kHZ
 F2=160; %sampling frequency of the original signals in the test set is 16 kHz
  
 output_dir='test';
 mkdir(output_dir);
 'load configuration file'
  load(namefile_configuration)

'simulate binaural signal'
order_conf=length(Az);

if order_conf==2 
    if nargin<3
    'eror syntax: namefile_configuration input_target input_noise'
    end
    
name=strcat(inputfile_target);
x1=load(name);x1=resample(x1,F1,F2);

name=strcat(inputfile_noise1);   
x2=load(name);x2=resample(x2,F1,F2);

N=min(length(x1),length(x2));

X=[reshape(x1(1:N),1,N);reshape(x2(1:N),1,N)];
end


if order_conf==3
    if nargin<4
    'eror syntax: namefile_configuration input_target input_noise1 input_noise2'
    end
    
name=strcat(inputfile_target);
x1=load(name);x1=resample(x1,F1,F2);

name=strcat(inputfile_noise1);   
x2=load(name);x2=resample(x2,F1,F2);

name=strcat(inputfile_noise2);   
x3=load(name);x3=resample(x3,F1,F2);

N=min([length(x1),length(x2),length(x3)]);

X=[reshape(x1(1:N),1,N);reshape(x2(1:N),1,N);reshape(x2(1:N),1,N)];
end

simulateSignal(X,Az,output_dir);

'feature extraction'

!extractFeatures test

'mask estimation'

estimateMask(output_dir,D,Region);

'resynthesis'

!resynthesize


%%%%%%%%%%END%%%%%%%%%%%