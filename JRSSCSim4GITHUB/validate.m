%% ValidateChanges

validationFile = 'output_files/testInputSerial.mat';
validationConfig = 'input_files/testInputSerial.m';
%validationConfig = 'input_files/testInputParpoolInner.m';

% Load to get 'saveOutput'
run(validationConfig);

% Run
ExpConsCheck(validationConfig);

% Check against known result.
% Changed order of operations to maximise potential difference.
r1 = load(validationFile, 'mle');
r2 = load(saveOutput, 'mle');
diff = sum(sum(abs(r1.mle - r2.mle)));
fprintf('diff=%g.\n', diff);
assert((diff < 1.e-6), 'Variation outside of acceptable limits.');
fprintf('Results look OK.');
