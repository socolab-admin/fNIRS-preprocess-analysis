% SYNTAX:
% params = create_params(save_flag)
%
%
% DESCRIPTION:
% This function initializes and returns a structured parameter object
% (`params`) used to configure an fNIRS preprocessing and analysis pipeline.
% The parameters define settings for trimming, channel quality assessment,
% motion artifact detection and correction, filtering, hemoglobin conversion,
% stimulus design, post-processing, and GLM HRF modelling.
%
% The function provides default values for all pipeline steps and allows
% optional saving of the parameter structure to a .mat file for reproducibility.
%
%
% INPUTS:
% save_flag: (optional, boolean)
%     If true, prompts the user to save the parameter structure to a .mat file.
%     If false or omitted, parameters are returned without saving.
%
%
% OUTPUTS:
% params: Struct `containing all pipeline configuration parameters.
%         Organized into the following fields:
%
% --- Paths ----
%   params.paths:
%       rawdir - Path to raw data directory.
%                Must contain group folders (e.g., G01, G02), each with
%                subject folders (e.g., S01_*), where each subject folder includes:
%                    - SNIRF file (*.snirf)
%                    - TRI file (*.tri)
%                    - probeInfo.mat
%       outdir - Path to output directory 
%
%   params.prefix:
%       group   - Prefix char used to identify group folders
%       subject - Prefix char used to identify subject folders
%
% --- Trimming ---
%   params.trim:
%       start_trigger - Trigger marking start of usable data; usually start
%                       of condition
%       end_trigger   - Trigger marking end of usable data; usually end
%                       of condition
%       buffer        - Time (sec) pre-stimulus onset (start trigger)
%                       Suggest adding a buffer when using GLM to properly align trange tPre [tPre tPost].
%                       Also helps account for filtering edge effects and HRF delay.
%
% --- Short Channels ---
%   params.short_channels:
%       map - Indices of short-separation channels (<~1–1.5 cm source-detector distance)
%             Fill in the correct values in square brackets []
%
% --- Homer version ---
%   params.nirs_use: Boolean flag indicating Homer2 (_Nirs) vs Homer3 processing.
%                    If true, Homer2 (_Nirs) functions will be used instead of Homer3 functions.
%                    Must be set to true if using motion correction method = "splineSG_custom".
%
%   ***run_preprocessing params below***
% --- Channel Quality Assessment ---
%   params.qc (quality control):
%   metrics    - Cell array of channel-quality metrics to compute: 
%                * "SCI" (scalp coupling index): quantifies optode–scalp coupling as the 
%                   zero-lag correlation between two wavelength signals 
%                   after band-pass filtering in the cardiac frequency band
%                   (e.g., 0.5–2.5 Hz) and normalization. 
%                   Values near 1 indicate strong shared cardiac pulsation 
%                   (good coupling), while values near 0 indicate poor signal quality 
%                   Citation: (Pollonini et al., 2013; Pollonini et al., 2016).

%                * "SNR" (signal-to-noise ratio): SNR quantifies channel signal quality 
%                   from raw light intensity in decibels (dB) as SNR(dB)=20log10(μ/σ), 
%                   where μ is the mean and σ is the standard deviation of the timeseries. 
%                   Higher values indicate more stable, higher-quality signals.
%                   Citation: (Yücel et al., 2021)

%                * "PP" (peak-power spectral ratio): PSP quantifies cardiac signal strength 
%                   by computing the peak spectral power of the cross-correlated, cardiac-band filtered, 
%                   and normalized wavelength signals within a sliding window (e.g., ~10 s). 
%                   Clean channels show a strong spectral peak at the cardiac frequency, 
%                   while motion-contaminated signals (possibly inflate SCI values) yield low or near-zero peak power. 
%                   Citation: (Pollonini et al., 2016). <still under
%                   development>

%                * "CHV" (cardiac-to-hemodynamic variability ratio): quantifies the relative variability 
%                   of the cardiac-band signal compared to the low-frequency hemodynamic-band signal 
%                   as the ratio of their standard deviations. 
%
%                You can list multiple metrics, e.g. ["SCI", "SNR"]. It is recommended to use multiple metrics to balance tradeoffs: 
%                for example, SCI reflects physiological signal quality, while SNR reflects overall signal stability (Yücel et al., 2021).
%
%   thresholds - Struct of thresholds for each metric: typical values in literature:
%                * SCI [DEFAULT 0.7] In adult populations under ideal experimental conditions, SCI thresholds are typically 
%                  set in the range of ~0.7–0.8 (Pollonini et al., 2014; Pollonini et al., 2016). However, SCI is sensitive 
%                  to factors such as motion, hair thickness, and cap fit, and lower thresholds (~0.5–0.6) are often used 
%                  for infants, children, and other special populations or more challenging datasets to avoid excessive channel rejection. 
%                  This flexibility reflects the need to balance data quality with data retention across different experimental contexts.

%                * SCI_band [DEFAULT [0.7, 1.5]] In young adult populations, ~50-90 bpm is common resting heart rate range. 
%                  Wider ranges (e.g., 0.5–2.5 Hz) are often recommended for more heterogeneous samples to better capture variability (Pollonini et al., 2013). 
%                  For infants and young children, higher heart rates justify broader ranges such as ~1.5–3 Hz (Richards, 1985; Fleming et al., 2011, Ostchega et al., 2011,). 
%                  In practice, users should tailor the SCI_band based on measured heart rate distributions when available, or adjust for specific
%                  populations and experimental conditions. Lastly to account for filtered edge effects, slightly widen band from
%                  proposed range.

%                * SNR [DEFAULT 20] ≤20 dB indicate higher relative noise and reduced signal reliability for acceptable fNIRS data quality (Yücel et al., 2021).

%                * PP [DEFAULT 0.1] spectral threshold to 0.10 accurately detects channels with strong cardiac pulsation that are free of motion artifacts (Pollonini et al., 2016).

%                * CHV [DEFAULT 1] In practice, ≥1 values often correspond to poor signal quality, reflecting noise or motion-dominated high-frequency 
%                  fluctuations and are typically associated with extremely low SCI values.

%                * motionClean [DEFAULT 0.6] ensures at least 60% of data is free from motion artifacts but adjust for infant and children populations. 

%                * badSubject [DEFAULT 0.5] A subject is marked as bad if the proportion of bad channels ≥ threshold.  
%                  Bad channels are defined as those failing selected QC metrics (e.g., SCI, SNR), exhibiting detector saturation or missing values (NaNs), or having insufficient  
%                  motion-free signal below threshold. In practice, subject exclusion cutoffs vary (~10–50% bad channels), reflecting tradeoffs between data quality and retention.
% 
%
% --- Motion artifact detection ---
%   params.motion.detect:
%       tMotion - [Default 1] Check for signal change indicative of a motion artifact over
%                 time range tMotion. Units of seconds.
%       tMask   - [Default 1] Mark data over +/- tMask seconds around the identified motion 
%                 artifact as a motion artifact. Units of seconds.
%       STDEV   - [Default 6] If the signal dod (delta optical density) for any given active channel changes by more
%                 that stdev_thresh * stdev(d) over the time interval tMotion, then
%                 this time point is marked as a motion artifact. The standard deviation is
%                 determined for each channel independently.
%       AMP     - [Default 5] If the signal dod for any given active channel changes by more
%                 that amp_thresh over the time interval tMotion, then this time point
%                 is marked as a motion artifact.
%
% --- Motion correction ---
%   params.motion.correct:
%       method - Selected motion correction method for
%        tPCA / rLOESS / wavelet / splineSG / splineSG_custom:
%       
%        * "tPCA" (Targetted PCA) - Identifies motion artifacts by flagging segments where any active channel
%                 exceeds STDEVthresh or AMPthresh; flagged segments (tIncAuto == 0) are
%                 then corrected using PCA filter.
%             nSV: Number of principal components to remove. If nSV < 1, components
%                  are removed until the specified fraction of variance is explained.
%                  Typical values (e.g., nSV = 0.97) remove ~97% of variance under the
%                  assumption that motion dominates global signal variance.
%             maxIter: Maximum number of recursive correction iterations. Common values
%                     (e.g., maxIter = 3) iteratively refine motion detection and
%                     correction until artifacts are minimized.
%             Citation: Yücel et al., 2014
%    
%        * "rLOESS" (robust locally weighted regression) - Applies robust locally weighted regression to smooth data and remove motion spikes.
%                   This nonparametric method fits local polynomial regressions with distance-based
%                   weighting and outlier resistance. This correction becomes computational expensive over many large-sample datasets.
%             span: Fraction of data used for local fitting; small values (e.g., span =\\\ 0.02) target
%                   fast, spike-like artifacts while preserving slower hemodynamic responses (<0.5 Hz).
%             Citation: Cleveland, 1979; Scholkmann et al., 2010
%
%       * "wavelet" (wavelet transform) - Perform a wavelet transformation of the dod data and computes the
%                 distribution of the wavelet coefficients. It sets the coefficient exceeding iqr times the interquartile range to zero, 
%                 because these are probably due to motion artifacts. This correction becomes computational expensive over many large-sample datasets.
%            iqr: 1.5 (why default)  parameter used to compute the statistics (iqr = 1.5 is 1.5 times the
%                 interquartile range and is usually used to detect outliers). Increasing it, it will delete fewer coefficients.
%                 If iqr<0 then this function is skipped. 
%            Citation: Molavi & Dumont, 2012; Cooper et al., 2012; Brigadoi et al., 2014
%
%        * "splineSG" (Spline + Savitzky-Golay) - Corrects motion artifacts by first
%                  estimating and removing baseline shifts using spline interpolation,
%                  followed by Savitzky-Golay filtering to smooth remaining spike-like
%                  artifacts. The Savitzky-Golay filter applies local polynomial (cubic)
%                  fitting over adjacent data points to preserve signal shape while
%                  reducing high-frequency noise.
%             p: Spline interpolation parameter controlling baseline correction
%                (default ≈ 0.99 as recommended in the literature). Lower values increase
%                correction strength. Set p < 0 to skip this method.
%             FrameSize_sec: Window length (in seconds) for the Savitzky-Golay filter,
%                            defining the size of the local fitting region.
%             Citation: Jahani et al., 2018 (Neurophotonics)
%
%        * "splineSG_custom" - Custom implementation of spline + Savitzky-Golay motion correction (not using Homer3 defaults).
%             p: Controls spline baseline correction strength.
%             FrameSize_sec: Defines smoothing window length (seconds).
%             k: Polynomial order used in Savitzky-Golay filtering.
%
%   params.motion.plot:
%       enable   - Enable/disable motion correction visualization;
%                  recommended only for debugging or small subsets of subjects.
%       channels - Channels to plot (column indices from SNIRF object);
%                  empty = automatically selects first 5 channels.
%
% --- Filtering ---
%   params.filter:
%       Applies bandpass filtering to fNIRS time series.
%       Default range (0.01–0.5 Hz) follows common recommendations for
%       removing slow drift and high-frequency noise (Yücel et al., 2021).
%       Narrower ranges (e.g., 0.01–0.1 Hz) are sometimes used for resting-state
%       or connectivity analyses to isolate low-frequency signals.
%       For GLM analyses, a higher low-pass cutoff (e.g., up to 0.5 Hz) is often
%       acceptable, as low-frequency drift is modeled separately via 3rd order polynomial terms.
%
%       use - Enable/disable filtering.
%       hpf - High-pass cutoff frequency (Hz); typically 0–0.01 to remove slow drift.
%       lpf - Low-pass cutoff frequency (Hz); typically 0.1–0.5 to remove physiological noise.
%
% --- Hemoglobin Concentration Conversion ---
%   params.hb:
%       ppf    - Partial pathlength factors (one per wavelength).
%                Commonly set to ~6 for adult head measurements assuming
%                homogeneous absorption (Duncan et al., 1996). Smaller values
%                reflect more localized absorption changes.
%                Increasingly, ppf = 1 is used (no scaling by pathlength or
%                source-detector distance), yielding units in Molar·mm (or cm);
%                this is a practical convention rather than a strict standard.
%       enable - Enable/disable hemoglobin concentration conversion.
%                Typically always enabled; disable only for debugging.
%       channels - Channels to visualize (indices from SNIRF object);
%                  empty = automatically selects first 5 channels.
%
% ***run_postprocessing params below***
% --- Postprocessing ---
%   params.post:
%       skip_bad   - Skip bad subjects if true.
%
%       use_short  - Include short-channel regressors for nuisance regression;
%                    improves removal of superficial physiological signals
%                    (e.g., scalp blood flow) and enhances recovery of cortical
%                    hemodynamics (Gagnon et al., 2011; Saager & Berger, 2008).
%
%       use_accel  - Include accelerometer regressors for nuisance regression;
%                    helps remove motion-related artifacts by modeling head
%                    movement contributions to the signal.
%
%       use_zscore - Apply within-channel z-scoring after nuisance regression
%                    to normalize signal amplitude across channels and subjects.
%       
% Stimulus Design
%   Defines condition timing used for epoching, block averaging, and GLM design matrix construction.
%   For multiple conditions, define each as cond(1), cond(2), ..., in numerical order.
%
%   params.post.stim:
%       cond(N):
%           name            - Condition name.
%           start_triggers  - Trigger codes marking condition onset.
%           end_triggers    - Trigger codes marking condition offset.
%           duration        - Duration of condition (sec) if no end trigger is used.
%
%   Example:
%       cond(1).name           = 'Prime';
%       cond(1).start_triggers = [1];
%       cond(1).end_triggers   = [2];
%       cond(1).duration       = 30;   % used when no end trigger is available
%
%       cond(2).name           = 'Video';
%       cond(2).start_triggers = [3];
%       cond(2).end_triggers   = [4];
%       cond(2).duration       = 90;   
% 
% Epoch
%   params.post.epoch:
%       use - Enable epoch extraction.
%       use_zscore - Z-score epoch window if true.
%       hrf_delay_sec - Hemodynamic response delay (sec) applied during epoching,
%                       shifting the onset/offset window to account for delayed
%                       Hb responses. Typical delays are ~4–5 s for HbO and
%                       ~5–7 s for HbR/BOLD signals (Huppert et al., 2006).
%
%
%
function params = create_params(save_flag)

% Create default params
params = struct();

%% ===============================
%% PATHS
%% ===============================
params.paths.rawdir = '';
params.paths.outdir = '';

params.prefix.group = 'G';
params.prefix.subject = 'S';

%% ===============================
%% TRIMMING
%% ===============================
params.trim.start_trigger = 3;   % replace with your own trigger code
params.trim.end_trigger   = 34;  % replace with your own trigger code
params.trim.buffer = 0;
% ===============================
%% SHORT CHANNEL INDICES
%% ===============================
params.short_channels.map = []; % add short channel indices (required for now)

%% ===============================
%% HOMER VERSION
%% ================================
params.nirs_use = false;  % true = Homer2 (_Nirs), false = Homer3

%% ===============================
%% CHANNEL QUALITY ASSESSMENT
%% ================================
params.qc.metrics = "SCI";
params.qc.thresholds.SCI = 0.70;
params.qc.thresholds.SCI_band = [0.7, 1.5];
params.qc.thresholds.SNR = 20.00;
params.qc.thresholds.CHV = 1.00;
params.qc.thresholds.PP = 0.10;
params.qc.thresholds.motionClean  = 0.60;
params.qc.thresholds.badSubject = 0.50;


%% ===============================
%% MOTION DETECTION
%% ================================
params.motion.detect.tMotion = 1;
params.motion.detect.tMask   = 1;
params.motion.detect.STDEV   = 6;
params.motion.detect.AMP     = 5;

%% ===============================
%% MOTION CORRECTION
%% ================================
params.motion.plot.enable = false; % leave this off as it will plot every subject's channel motion corrected signals (be careful to plot subject data you want to see)         % turn on/off
params.motion.plot.channels = [];   % channel idicies you want to plot; if leave empty then would plot first 5 from starting 0 indx
params.motion.correct.method = "tPCA";   % "rLOESS", "wavelet", "splineSG", "splineSG_custom", "none"

% ---- tPCA ----
params.motion.correct.tPCA.nSV = 0.9;
params.motion.correct.tPCA.maxIter = 5;

% ---- rLOESS ----
params.motion.correct.rLOESS.span = 0.02;

% ---- Wavelet ----
params.motion.correct.wavelet.iqr = 1.5;

% ---- Homer3 Spline + SG ----
params.motion.correct.splineSG.p = 0.99;
params.motion.correct.splineSG.FrameSize_sec = 10;

% ---- Custom (Homer2-style) ----
params.motion.correct.splineSG_custom.p = 0.99;
params.motion.correct.splineSG_custom.FrameSize_sec = 2;
params.motion.correct.splineSG_custom.K = 3;

%% ===============================
%% HEMOGLOBIN CONCENTRATION
%% ================================
params.hb.plot.enable = false;
params.hb.plot.channels = []; 
params.hb.ppf = [6 6];

%% ===============================
%% FILTER
%% ================================
params.filter.use = true;
params.filter.hpf = 0.01;
params.filter.lpf = 0.5;

%% ===============================
%% POST-PROCESSING (NUISANCE, ZSCORE)
%% ===============================
params.post.skip_bad = false;

% Stim Design Matrix
% if you want to add more conditions just drop down params.post.stim.cond(N).name/triggers/duration 
% if no params for this then this is skipped and no stims is built for data
params.post.stim.cond(1).name = 'Cond1';  % replace with actual condition name
params.post.stim.cond(1).start_triggers = [];   % add trigger code onsets
params.post.stim.cond(1).end_triggers = []; % add trigger code onsets 
% make sure to have start and end triggers balanced in number of trigger
% codes etc.
params.post.stim.cond(1).duration = 30; % replace with actual stimulus duration 

% Nuisance regressors
params.post.use_short = false;     % use short channels
params.post.use_accel = false;    % use accelerometer

% Z-scoring
params.post.use_zscore = false;

% Epochs
params.post.epoch.use = false;
params.post.epoch.use_zscore = false;
params.post.epoch.hrf_delay_sec = 0; % 4.5 for specific canonical HbO peak and 6 for canonical HbR peak; 6 conservative peak for both concentrations

% If saves parameters
if nargin > 0 && save_flag == true

    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    default_name = ['params_' timestamp '.mat'];

    [file, path] = uiputfile('*.mat', 'Save parameter file as', default_name);

    if isequal(file,0)
        disp('User canceled saving params.');
    else
        save(fullfile(path, file), 'params');
        fprintf('Params saved to: %s\n', fullfile(path, file));
    end
end

end