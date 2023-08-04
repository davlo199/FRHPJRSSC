validationFile = 'output_files/testScaleBig.mat';
validationConfig = 'input_files/testScaleBig.m';
% Load to get 'saveOutput'
run(validationConfig);
% Run
ExpConsCheck(validationConfig);
