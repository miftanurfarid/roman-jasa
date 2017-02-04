%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%TRAINING for a configuration with two or three sources %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function run_training(Az,namefile_configuration)

% Input: Az=[azimuth_target,azimuth_interference]
% Useful output files: namefile_configuration to store azimuth information and decision rules.
 
% %%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%TARGET ITD%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%

mkdir ('target');

x = randn(1,441*40)*3000; %400 ms of white noise simulated at target location
simulateSignal(x,Az(1),'target');

if ispc == 1
    !targetITD target
end

if isunix == 1
    !./targetITD target
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%BINAURAL SIMULATION AND FEATURE EXTRACTION%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

order_conf = length(Az);

if order_conf < 2 || order_conf > 3
    disp('Configuration should have two or three source locations');
end

if order_conf == 2
    F1 = 441; %sampling frequency of the HRTF database is 44.1 kHZ
    F2 = 160; %sampling frequency of the original training signals is 16 kHz
    
    mkdir ('train');
    output_dir = 'train';

    for ind_T = 0:4
        for ind_N = 5:9
            fprintf('ind_T = %d\n',ind_T)
            fprintf('ind_N = %d\n',ind_N)

            name = strcat('TRAINSET\','T',int2str(ind_T));
            x1   = load(name);
            x1   = resample(x1,F1,F2);

            name = strcat('TRAINSET\','T',int2str(ind_N));   
            x2   = load(name);
            x2   = resample(x2,F1,F2);

            N = min(length(x1),length(x2));

            X = [x1(1:N);x2(1:N)];

            simulateSignal(X,Az,output_dir);

            if ispc == 1
                !extractFeatures train 
            end

            if isunix == 1
                !./extractFeatures train 
            end
        end
    end
end
    

if order_conf == 3
    F1 = 441; %sampling frequency of the HRTF database is 44.1 kHZ
    F2 = 160; %sampling frequency of the original training signals is 16 kHz
    
    mkdir train;
    output_dir='train';
    
    for ind_T = 0:3
        for ind_N1 = 4:6
            for ind_N2 = 7:9
                
                fprintf('ind_T %d\n',ind_T);
                fprintf('ind_N1 %d\n',ind_N1);
                fprintf('ind_N2 %d\n',ind_N2);
                
                name = strcat('TRAINSET\','T',int2str(ind_T));
                x1   = load(name);
                x1   = resample(x1,F1,F2);
                
                name = strcat('TRAINSET\','T',int2str(ind_N1));
                x2   = load(name);
                x2   = resample(x2,F1,F2);
                
                name = strcat('TRAINSET\','T',int2str(ind_N2));   
                x3   = load(name);
                x3   = resample(x3,F1,F2);
                
                N    = min([length(x1),length(x2),length(x3)]);
                
                X    = [x1(1:N);x2(1:N);x3(1:N)];
                
                simulateSignal(X,Az,output_dir);

                if ispc == 1
                    !extractFeatures train 
                end

                if isunix == 1
                    !./extractFeatures train 
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%% CLASSIFICATION%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

training('ITD', 'IID', 'RATIO', Az, namefile_configuration)
 
%%%%%END TRAINING%%%%%%%

end