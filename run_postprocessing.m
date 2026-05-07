% SYNTAX:
% data_postprocessed = run_postprocessing(params)
%
%
% DESCRIPTION:
% This function performs postprocessing of preprocessed fNIRS data (.mat files)
% organized by group and subject folders. It includes subject-level filtering
% (e.g., skipping bad subjects), construction of stimulus design matrices,
% nuisance regression (short-channel and/or accelerometer), optional z-scoring,
% and epoch extraction. Outputs are saved per subject and aggregated across subjects.
%
%
% INPUTS:
% params: Structure containing all postprocessing settings and paths.
%         Parameters are defined in create_params().
%
%
% INPUT FILE REQUIREMENTS:
% Each subject folder must contain:
%     - Preprocessed file: [subjectID '_preprocessed.mat']
%       (output from run_preprocessing)
%
%
% OUTPUTS:
% data_postprocessed - Cell array containing all postprocessed subject data.
%                      Each cell corresponds to a subject and includes:
%                          - preprocessed.mat fields
%                          - HbO, HbR (postprocessed signals)
%                          - stim (condition design matrix)
%                          - subjectID
%                          - post.epochs (if epoching enabled)
%
% For each subject, a postprocessed .mat file is also saved:
%
% outfile = [subjectID '_postprocessed.mat']
%
% containing:
%     data_postprocessed - Struct for that subject only
%
% Additionally, per-epoch CSV files are saved (if enabled):
%     subjectID_condition_postprocessed_EpochXX_HbO.csv
%     subjectID_condition_postprocessed_EpochXX_HbR.csv
%
%
% OUTPUT DIRECTORY STRUCTURE:
%     outdir / groupID / subjectID / subjectID_postprocessed.mat
%     outdir / groupID / subjectID / subjectID_*_EpochXX_*.csv
%
%
% POSTPROCESSING STEPS:
% 1. Load preprocessed subject data
% 2. Optionally skip subjects flagged as bad (QC-based)
% 3. Identify bad channels and short-separation channels
% 4. Construct stimulus design matrix from trigger information
% 5. Build nuisance regressors:
%       - Short-channel signals (if enabled)
%       - Accelerometer signals (if available and enabled)
% 6. Perform nuisance regression on long channels
% 7. Optionally apply within-channel z-scoring
% 8. Reconstruct full channel matrix (short + long channels)
% 9. Extract epochs per condition with optional HRF delay shift
% 10. Save postprocessed outputs and epoch files
%
%
% NOTES:
%     - Requires outputs from run_preprocessing (HbO, HbR, QC, triggers, etc.).
%     - Short-channel regression removes superficial physiological noise.
%     - Accelerometer regressors help account for motion-related artifacts.
%     - Epoch extraction applies an HRF delay to align neural responses.
%     - Z-scoring is applied after regression to normalize signal amplitudes.
%     - data_postprocessed (function output) aggregates all subjects, while
%       each saved .mat file contains only that subject’s data.
function data_postprocessed = run_postprocessing(params)

data_postprocessed = {};

%% ===============================
%% PATHS
%% ===============================

outdir = params.paths.outdir;

groupprefix   = params.prefix.group;
subjectprefix = params.prefix.subject;

groupdirs = dir(fullfile(outdir, [groupprefix '*']));
groupdirs = groupdirs([groupdirs.isdir]);

for i = 1:length(groupdirs)

    groupname = groupdirs(i).name;
    groupfolder = fullfile(outdir, groupname);

    subjectdirs = dir(fullfile(groupfolder, [subjectprefix '*']));
    subjectdirs = subjectdirs([subjectdirs.isdir]);

    for j = 1:length(subjectdirs)

        subjectname = subjectdirs(j).name;
        subjectfolder = fullfile(groupfolder, subjectname);

        subjTok = regexp(subjectname, '^(.*?)_', 'tokens', 'once');
        if isempty(subjTok)
            continue
        end
        subjectID = subjTok{1};

        scanfile = fullfile(subjectfolder, [subjectID '_preprocessed.mat']);

        if ~isfile(scanfile)
            fprintf('Missing file, skipping: %s\n', scanfile);
            continue
        end
        
        data_preprocessed = load(scanfile);

        % ===============================
        % Skip bad subjects
        % ===============================
        is_bad_subject = any(table2array(data_preprocessed.quality_report("BadSub", :)));

        if params.post.skip_bad && is_bad_subject
            fprintf('Skipping bad subject %s\n', subjectID);
            continue
        end

        fprintf('Processing subject %s\n', subjectID);

        % ===============================
        % Extract data
        % ===============================
        HbO_full = data_preprocessed.HbO;
        HbR_full = data_preprocessed.HbR;

        nChan = size(HbO_full,2);

        % Channel masks
        idx_short = data_preprocessed.idx_short(:)';
        is_short = false(1,nChan);
        is_short(idx_short) = true;

        %is_bad_chan = false(1, size(data_preprocessed.quality_report,2));
        is_bad_chan = table2array(data_preprocessed.quality_report("BadChan", :)) == 1;

        %disp(is_bad_chan)
        fprintf('Detected %d bad channels\n', sum(is_bad_chan));

        %% ===============================
        %% Stim - Condition Design Matrix
        %% ===============================
        
        stim = StimClass();
        
        if isfield(params.post, 'stim') && isfield(params.post.stim, 'cond')
        
            t = data_preprocessed.t;
            fs = data_preprocessed.fs;
            tri_samples = data_preprocessed.tri_raw{2};
            tri_codes   = data_preprocessed.tri_raw{3};
        
            nCond = length(params.post.stim.cond);
        
            for c = 1:nCond
        
                cond = params.post.stim.cond(c);
        
                % ---- find trigger indices ----
                idx = find(ismember(tri_codes, cond.start_triggers));

                if isempty(idx)
                    warning('No triggers found for %s', cond.name);
                    continue
                end
        
                % ---- convert to time ----
                onset_samples = tri_samples(idx);
                start_idx = find(tri_codes == params.trim.start_trigger, 1, 'first');
                buffer = round(params.trim.buffer*fs); 
                %disp(buffer)
                firstSample = tri_samples(start_idx) - buffer;
                %start_idx = find(ismember(tri_codes, params.trim.start_trigger));
                onset_samples = onset_samples - firstSample + 1;
                %disp(onset_samples);
                cond_onsets = t(onset_samples);
                %disp(cond_onsets)

                % ---- duration + amplitude ----
                n = length(cond_onsets);
        
                durations  = cond.duration * ones(n,1);
                amplitudes = ones(n,1);
        
                % ---- build stim ----
                stim(c).SetName(cond.name);
                stim(c).SetData([cond_onsets(:), durations, amplitudes]);
        
                % ---- debug ----
                fprintf('Condition %s: %d events\n', cond.name, n);
                end
            
            else
                warning('No stim parameters found — skipping stim creation');
        end
            
        %disp(stim)
        data_preprocessed.stim = stim;

        %% ===============================
        %% Process HbO / HbR
        %% ===============================
        HbO_proc = process_Hb_simple(params, data_preprocessed, HbO_full, is_bad_chan, is_short);
        HbR_proc = process_Hb_simple(params, data_preprocessed, HbR_full, is_bad_chan, is_short);

        %% ===============================
        %% EPOCH 
        %% ===============================
        if params.post.epoch.use
            data_preprocessed = epochs(params, data_preprocessed, subjectfolder, subjectID, HbO_proc, "HbO");
            data_preprocessed = epochs(params, data_preprocessed, subjectfolder, subjectID, HbR_proc, "HbR");
        end

        % ===============================
        % SAVE
        % ===============================
        data_preprocessed.subjectID = subjectID;
        data_preprocessed.HbO = HbO_proc;
        data_preprocessed.HbR = HbR_proc;
        
        data_postprocessed{i,j} = data_preprocessed;

        outfile = fullfile(subjectfolder, [subjectID '_postprocessed.mat']);
        save(outfile, 'data_postprocessed');

    end
end
end

function Hb_out = process_Hb_simple(params, data_preprocessed, Hb_full, is_bad_chan, is_short)

%% ===============================
%% Mask bad channels
%% ===============================
Hb = Hb_full;
Hb(:, is_bad_chan) = NaN;

%% ===============================
%% Split channels
%% ===============================
Hb_short = Hb(:, is_short);
Hb_long  = Hb(:, ~is_short);

%% ===============================
%% Build nuisance regressors
%% ===============================
Hb_nuisance = build_nuisance(params, data_preprocessed, Hb_short, size(Hb,1));

% Remove NaNs
Hb_nuisance(:, any(isnan(Hb_nuisance),1)) = [];

%% ===============================
%% Nuisance regression
%% ===============================
Hb_long_reg = Hb_long;

for ch = 1:size(Hb_long,2)
    if ~isnan(Hb_long(1,ch))
        [~,~,r] = regress(Hb_long(:,ch), Hb_nuisance);
        Hb_long_reg(:,ch) = r;
    end
end

%% ===============================
%% Z-score 
%% ===============================
if params.post.use_zscore
    Hb_long_reg = zscore(Hb_long_reg);
    Hb_short    = zscore(Hb_short);
end

%% ===============================
%% Reconstruct
%% ===============================
Hb_out = NaN(size(Hb));

Hb_out(:, ~is_short) = Hb_long_reg;
Hb_out(:,  is_short) = Hb_short;

end

function Hb_nuisance = build_nuisance(params, data_preprocessed, Hb_short, T)

Hb_nuisance = ones(T,1);  % intercept

if params.post.use_short
    Hb_nuisance = [Hb_nuisance, Hb_short];
end

if params.post.use_accel && isfield(data_preprocessed,'accelerometer') && ~isempty(data_preprocessed.accelerometer)
    Hb_nuisance = [Hb_nuisance, data_preprocessed.accelerometer];
end

end

function data_preprocessed = epochs(params, data_preprocessed, subjectfolder, subjectID, Hb, HbType)

stim = data_preprocessed.stim;
t = data_preprocessed.t;
fs = data_preprocessed.fs;
tri_samples = data_preprocessed.tri_raw{2};
tri_codes   = data_preprocessed.tri_raw{3};

HRF_shift = round(params.post.epoch.hrf_delay_sec * fs);
%disp(HRF_shift);

chanSrcDet = data_preprocessed.probeInfo.probes.index_c;
chanLabels = arrayfun(@(i) ...
    sprintf('S%d-D%d', chanSrcDet(i,1), chanSrcDet(i,2)), ...
    1:size(chanSrcDet,1), 'UniformOutput', false);

labels = ['time', chanLabels(:)'];

for c = 1:length(stim)

    condName = params.post.stim.cond(c).name;

    % ===============================
    % ONSETS
    % ===============================
    idx = ismember(tri_codes, params.post.stim.cond(c).start_triggers);
    onsets_samples = tri_samples(idx);
    start_idx = find(tri_codes == params.trim.start_trigger, 1, 'first');
    buffer = round(params.trim.buffer*fs); 
    firstSample = tri_samples(start_idx) - buffer;
    onsets_samples = onsets_samples - firstSample + 1;
    %disp(onsets_samples);

    % ===============================
    % OFFSETS
    % ===============================
    idx = ismember(tri_codes, params.post.stim.cond(c).end_triggers);
    offsets_samples = tri_samples(idx);
    start_idx = find(tri_codes == params.trim.start_trigger, 1, 'first');
    buffer = round(params.trim.buffer*fs); 
    firstSample = tri_samples(start_idx) - buffer;
    %disp(firstSample);
    offsets_samples = offsets_samples - firstSample + 1;
    %disp(offsets_samples);

    % --- Number of epochs ---
    nEpochs = length(onsets_samples);

    %offsets_sec = t(offsets_samples);
    %print("print offsets sec")
    %disp(offsets_sec)
    %onsets_idx = round(onsets_sec * fs);
    %disp(onsets_idx)
    %offset_idx = round(offsets_sec * fs);
    %disp(offsets_idx)

    for k = 1:nEpochs

        s = onsets_samples(k) + HRF_shift;
        %disp(s)
        e = offsets_samples(k) + HRF_shift;
        %disp(e)

        if s < 1 || e > size(Hb,1)
            continue
        end

        block = Hb(s:e,:);
        %% ===============================
        %% Z-score Epochs
        %% ===============================
        if params.post.epoch.use_zscore
            block = zscore(block);
        end
        T = array2table([data_preprocessed.t(s:e), block], 'VariableNames', labels);
        
        % ===============================
        % SAVE FILE
        % ===============================
        outfile = fullfile(subjectfolder, ...
            sprintf('%s_%s_postprocessed_Epoch%02d_%s.csv', subjectID, condName, k, HbType));

        writetable(T, outfile);

        % ===============================
        % SAVE TO STRUCT
        % ===============================
        data_preprocessed.post.epochs.(HbType).(condName){k} = T;

       % save / store
    end
end

end

% ===============================
% HELPER FUNCTIONS 
% ===============================

function thresh = get_thresh(params, name, default_val)

if isfield(params.qc.thresholds, name)
    thresh = params.qc.thresholds.(name);
else
    thresh = default_val;
end

end

function val = get_motion_param(params, method, field, default_val)

if isfield(params.motion.correct, method)
    if isfield(params.motion.correct.(method), field)
        val = params.motion.correct.(method).(field);
        return
    end
end

val = default_val;

end

