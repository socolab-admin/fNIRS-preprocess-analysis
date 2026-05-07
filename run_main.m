clear; clc;

% Step 1: create params
params = create_params();

% Step 2 (Optional): set required values if haven't in params (uncomment
% code to search for rawdir and outdir
params.paths.rawdir = uigetdir('', 'Select RAW data folder');
params.paths.outdir = uigetdir('', 'Select OUTPUT folder');

% Step 3: save params log file used for methods clarity
outfile = fullfile(params.paths.outdir, 'params_summary.xlsx');
save_params_to_excel(params, outfile);

% Step 4: run preprocessing
run_preprocessing(params);
    
% Step 5: run postprocessing
data = run_postprocessing(params);

% Step 6: run glm
run_glm(data);