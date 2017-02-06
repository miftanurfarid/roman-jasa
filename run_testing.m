%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%TESTING EXAMPLE - TWO-SOURCE CONFIGURATION%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function run_testing(namefile_configuration, inputfile_target, inputfile_noise1,inputfile_noise2)

F1 = 441; %sampling frequency of the HRTF database is 44.1 kHZ
  
output_dir='test';
mkdir(output_dir);
disp('load configuration file');
load(namefile_configuration);

disp('simulate binaural signal');
order_conf = length(Az);

if order_conf == 2
    if nargin < 3
    disp('error syntax: namefile_configuration input_target input_noise')
    end
    
    name = strcat(inputfile_target);
    [x1,fs] = audioread(name);
    x1 = resample(x1,F1,fs);
    
    name = strcat(inputfile_noise1);
    [x2,fs] = audioread(name);
    x2 = resample(x2,F1,fs);
    
    N = min(length(x1),length(x2));
    
    X = [reshape(x1(1:N),1,N);reshape(x2(1:N),1,N)];
end


if order_conf == 3
    if nargin < 4
        disp('error syntax: namefile_configuration input_target input_noise1 input_noise2')
    end
    
    name = strcat(inputfile_target);
    [x1,fs] = audioread(name);
    x1 = resample(x1,F1,fs);
    
    name = strcat(inputfile_noise1);
    [x2,fs] = audioread(name);
    x2 = resample(x2,F1,fs);

    name = strcat(inputfile_noise1);
    [x3,fs] = audioread(name);
    x3 = resample(x3,F1,fs);
    
    N = min([length(x1),length(x2),length(x3)]);
    
    X = [reshape(x1(1:N),1,N);reshape(x2(1:N),1,N);reshape(x2(1:N),1,N)];
end

simulateSignal(X,Az,output_dir);

disp('feature extraction');

if ispc == 1
    !extractFeatures test
end

if isunix == 1
    !./extractFeatures test
end

disp('mask estimation');

estimateMask(output_dir,D,Region);

disp('resynthesis');

if ispc == 1
    !resynthesize
end

if isunix == 1
    !./resynthesize
end

end
%%%%%%%%%%END%%%%%%%%%%%