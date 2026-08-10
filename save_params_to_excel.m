% SYNTAX:
% save_params_to_excel(params, outfile)
%
%
% DESCRIPTION:
% Saves preprocessing and postprocessing parameters to an Excel file,
% organized by pipeline steps for documentation and reproducibility.
%
%
% INPUTS:
% params  - Structure containing all pipeline parameters (see create_params()).
% outfile - Output Excel file path.
%
%
% OUTPUTS:
% Excel file containing ordered parameters grouped by step
% (trimming, QC, motion, filtering, Hb conversion, postprocessing, etc.).
%
%
% NOTES:
%     - Only selected QC metrics and chosen motion method parameters are saved.
%     - Useful for tracking and sharing pipeline configurations.
function save_params_to_excel(params, outfile)

rows = {};

% helper
function add(order, step, param, value)
    rows(end+1,:) = {order, step, param, value};
end

%% ===============================
% 1. TRIMMING
% ===============================
add(1, 'Trimming', 'Start Trigger', params.trim.start_trigger);
add(1, 'Trimming', 'End Trigger', params.trim.end_trigger);
add(1, 'Trimming', 'Buffer (sec)', params.trim.buffer);

%% ===============================
% 2. CHANNEL QUALITY (SELECTED ONLY)
% ===============================
selected = params.qc.metrics;

add(2, 'Channel QC', 'Selected Metrics', strjoin(selected, ', '));

for i = 1:length(selected)

    metric = selected(i);

    switch metric

        case "SCI"
            add(2, 'Channel QC', 'SCI Threshold', params.qc.thresholds.SCI);

        case "SNR"
            add(2, 'Channel QC', 'SNR Threshold', params.qc.thresholds.SNR);

        case "QT"
            add(2, 'Channel QC', 'QT Threshold', params.qc.thresholds.QT);

        case "PP"
            add(2, 'Channel QC', 'PSP Threshold', params.qc.thresholds.PSP);

    end
end

add(2, 'Channel QC', 'Bad Channel Prop Threshold', params.qc.thresholds.badSubject);

%% ===============================
% 3. MOTION DETECTION
% ===============================
add(3, 'Motion Detection', 'tMotion', params.motion.detect.tMotion);
add(3, 'Motion Detection', 'tMask', params.motion.detect.tMask);
add(3, 'Motion Detection', 'STDEV Threshold', params.motion.detect.STDEV);
add(3, 'Motion Detection', 'AMP Threshold', params.motion.detect.AMP);

%% ===============================
% 4. MOTION CORRECTION (SELECTED ONLY)
% ===============================
method = params.motion.correct.method;

add(4, 'Motion Correction', 'Method', method);

switch method

    case "tddr"
        add(4, 'Motion Correction', 'Note', 'tddr requires no parameters');

    case "tPCA"
        add(4, 'Motion Correction', 'nSV', params.motion.correct.tPCA.nSV);
        add(4, 'Motion Correction', 'maxIter', params.motion.correct.tPCA.maxIter);

    case "rhloess"
        add(4, 'Motion Correction', 'span', params.motion.correct.rLOESS.span);

    case "wavelet"
        add(4, 'Motion Correction', 'iqr', params.motion.correct.wavelet.iqr);

    case "splineSG"
        add(4, 'Motion Correction', 'p', params.motion.correct.splineSG.p);
        add(4, 'Motion Correction', 'FrameSize_sec', params.motion.correct.splineSG.FrameSize_sec);

    case "splineSG_custom"
        add(4, 'Motion Correction', 'p', params.motion.correct.splineSG_custom.p);
        add(4, 'Motion Correction', 'FrameSize_sec', params.motion.correct.splineSG_custom.FrameSize_sec);
        add(4, 'Motion Correction', 'K', params.motion.correct.splineSG_custom.K);

    case "none"
        add(4, 'Motion Correction', 'Note', 'No correction applied');

end

%% ===============================
% 5. FILTERING
% ===============================
add(5, 'Filtering', 'Enabled', params.filter.use);
add(5, 'Filtering', 'Low cutoff (Hz)', params.filter.hpf);
add(5, 'Filtering', 'High cutoff (Hz)', params.filter.lpf);

%% ===============================
% 6. HB CONVERSION
% ===============================
add(6, 'Hb Conversion', 'PPF', mat2str(params.hb.ppf));

%% ===============================
% 7. POST-PROCESSING
% ===============================
add(7, 'Post', 'Regress Short Channels', params.post.use_short);
add(7, 'Post', 'Regress Accelerometer', params.post.use_accel);
add(7, 'Post', 'Z-score', params.post.use_zscore);

%% ===============================
% 8. STIM DESIGN
% ===============================
if isfield(params.post, 'stim')
    for c = 1:length(params.post.stim.cond)
        cond = params.post.stim.cond(c);

        add(8, 'Stim', sprintf('Cond%d Name', c), cond.name);
        add(8, 'Stim', sprintf('Cond%d Start Triggers', c), mat2str(cond.start_triggers));
        add(8, 'Stim', sprintf('Cond%d End Triggers', c), mat2str(cond.end_triggers));
        add(8, 'Stim', sprintf('Cond%d Duration', c), cond.duration);
    end
end

%% ===============================
% 9. EPOCH
% ===============================
add(9, 'Epoch', 'Enabled', params.post.epoch.use);
add(9, 'Epoch', 'Z-score', params.post.epoch.use_zscore);
add(9, 'Epoch', 'HbO-HRF Delay (sec)', params.post.epoch.hbo_hrf_delay_sec);
add(9, 'Epoch', 'HbR-HRF Delay (sec)', params.post.epoch.hbr_hrf_delay_sec);

%% ===============================
% BUILD TABLE
% ===============================
T = cell2table(rows, ...
    'VariableNames', {'Order', 'Step', 'Parameter', 'Value'});

%% ===============================
% SAVE
% ===============================
writetable(T, outfile);

fprintf('Params saved to Excel: %s\n', outfile);

end
