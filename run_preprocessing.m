% SYNTAX:
% run_preprocessing(params)
%
%
% DESCRIPTION:
% This function performs end-to-end preprocessing of fNIRS SNIRF data organized
% by group and subject folders. It includes trimming using trigger (.tri) files,
% quality control (QC) assessment, motion artifact detection and correction,
% conversion to hemoglobin concentration (HbO, HbR, HbT), optional filtering,
% and data visualization. Processed outputs are saved per subject.
%
%
% INPUTS:
% params: Structure containing all preprocessing settings and paths.
%         Parameters are defined in create_params().
%
% INPUT FILE REQUIREMENTS:
% Each subject folder must contain:
%     - SNIRF file (*.snirf): raw fNIRS data
%     - TRI file (*.tri): trigger timestamps and codes for trimming
%     - probeInfo.mat: source-detector mapping information
%
%
% OUTPUTS:
% For each subject, a preprocessed .mat file is saved:
%
% outfile = [subjectID '_preprocessed.mat']
%
% containing:
%     HbO            - Oxygenated hemoglobin time series
%     HbR            - Deoxygenated hemoglobin time series
%     HbT            - Total hemoglobin time series
%     dhbFilt        - Filtered hemoglobin data structure
%     tri_raw        - Raw trigger information
%     idx_short      - Short-separation channel indices
%     accelerometer  - Interpolated and filtered accelerometer data (if available)
%     quality_report - Channel-wise QC metrics table
%     fs             - Sampling rate (Hz)
%     t              - Time vector (sec)
%     SD             - Source-detector structure
%     probeInfo      - Probe geometry and channel mapping
%
% Additionally, a group-level QC summary is saved:
%
% group_qc_summary.mat containing:
%     group_qc.channel - Table of channel-wise bad rates across subjects
%     group_qc.metric  - Table of aggregated QC metric means
%     group_qc.subject - Table of subject-level pass/fail flags
%
% OUTPUT DIRECTORY STRUCTURE:
%     outdir / groupID / subjectID / subjectID_preprocessed.mat
%     outdir / group_qc_summary.mat
%
%
% PREPROCESSING STEPS:
% 1. Load SNIRF and TRI files
% 2. Trim data using trigger codes
% 3. Detect short channels and accelerometer signals
% 4. Compute QC metrics:
%       - Detector saturation
%       - Motion-free proportion
%       - SNR (dB)
%       - Scalp Coupling Index (SCI)
%       - SD ratio (cardiac vs hemodynamic)
%       - Bad channel and subject classification
% 5. Convert raw light intensity to delta optical density.
% 6. Detect motion artifacts 
% 6. Apply motion correction (tPCA, wavelet, splineSG, rLOESS, etc.)
% 7. Convert delta optical density → delta hemoglobin concentration
% 8. Apply optional bandpass filtering
% 9. Generate optional QC and signal plots
% 10. Save preprocessed outputs
%
%
%
% NOTES:
%     - Pipeline supports both Homer2-style and Homer3-style processing.
%     - Motion correction method is user-selectable via params.motion.correct.method.
%     - Accelerometer data (if present) is interpolated to fNIRS time base and can
%       be used for QC or regression.
%     - QC summary is generated at the group level after processing all subjects.
function run_preprocessing(params)

%% ===============================
%% PATHS
%% ===============================

rawdir = params.paths.rawdir;
outdir = params.paths.outdir;

groupprefix   = params.prefix.group;
subjectprefix = params.prefix.subject;

% -------------------------
% Find group folders
% -------------------------
groupdirs = dir(fullfile(rawdir, [groupprefix '*']));
groupdirs = groupdirs([groupdirs.isdir]);

% ---- QC data summary acculumators ---
metrics = ["DetectorSat", "MotionClean", "SNR", "SCI", "PSP", "QT", "BadChan", "BadSub"];
subject_fail_count = 0;
subject_total_count = 0;
all_quality_reports = {};   % store each subject's QC table
all_subject_ids = {};

for i = 1:length(groupdirs)

    groupname = groupdirs(i).name;
    groupfolder = fullfile(rawdir, groupname);

    % Extract clean group ID (e.g., G01)
    grpTok = regexp(groupname, '^(G\d+)', 'tokens', 'once');
    if isempty(grpTok)
        continue
    end
    groupID = grpTok{1};

    % -------------------------
    % Find subject folders
    % -------------------------
    subjectdirs = dir(fullfile(groupfolder, [subjectprefix '*']));
    subjectdirs = subjectdirs([subjectdirs.isdir]);

    for j = 1:length(subjectdirs)

        subjectname = subjectdirs(j).name;
        subjectfolder = fullfile(groupfolder, subjectname);

        % Extract clean subject ID (e.g., S01)
        subjTok = regexp(subjectname, '^(.*?)_', 'tokens', 'once');
        if isempty(subjTok)
            continue
        end
        subjectID = subjTok{1};

        % -------------------------
        % Find SNIRF files
        % -------------------------
        snirf_files = dir(fullfile(subjectfolder, '*.snirf'));

        if isempty(snirf_files)
            warning('No SNIRF file in %s', subjectfolder);
            continue
        end

        % -------------------------
        % Check for TRI file
        % -------------------------
        tri_files = dir(fullfile(subjectfolder, '*.tri'));
        has_tri = ~isempty(tri_files);

        % -------------------------
        % Print subject info
        % -------------------------
        fprintf('\nGroup:   %s\n', groupID);
        fprintf('Subject: %s\n', subjectID);

        if ~has_tri
            fprintf('SKIPPING: Missing .tri file (no pre-processing)\n');
            continue   
        end

        % -------------------------
        % Load probeInfo (channel mapping)
        % -------------------------
        probeInfoFile = dir(fullfile(subjectfolder, '*probeInfo.mat'));
        
        if isempty(probeInfoFile)
            warning('Missing probeInfo.mat in %s', subjectfolder);
            continue
        end
        
        tmp = load(fullfile(subjectfolder, probeInfoFile(1).name));
        probeInfo = tmp.probeInfo;
        
        chanSrcDet = probeInfo.probes.index_c;   % [Nchannels x 2]

        %% 
        % -------------------------
        % Process SNIRF files
        % -------------------------
        for k = 1:length(snirf_files)

            snirfname = snirf_files(k).name;
            snirfpath = fullfile(subjectfolder, snirfname);

            % Create output directory
            outpath = fullfile(outdir, groupname, subjectname);
            if ~exist(outpath, 'dir')
                mkdir(outpath)
            end

            fprintf('SNIRF:   %s\n', snirfpath);

            
            % ---- Start preprocessing ----
            snirf = SnirfLoad(snirfpath);
            
            %% ===============================
            %% TRIMMING
            %% ===============================
            %snirf_trimmed = snirf.copy;
            snirf_trimmed.data = snirf.data.copy;
            snirf_trimmed.aux = snirf.aux.copy;
            snirf_trimmed.probe = snirf.probe.copy;
            snirf_trimmed.stim = snirf.stim.copy;
            snirf_trimmed.metaDataTags = snirf.metaDataTags.copy;

            % ---- Load .tri file ----
            tri_path = fullfile(subjectfolder, tri_files(1).name);

            fid = fopen(tri_path, 'r');
            tri_raw = textscan(fid, '%s %f %f', 'Delimiter', ';');
            fclose(fid);

            tri_time_str = tri_raw{1};   %#ok<NASGU> % timestamp strings
            tri_samples  = tri_raw{2};   % sample indices
            tri_codes    = tri_raw{3};   % trigger codes

            %disp(tri_samples);
            %disp(tri_codes);


            % ---- Sampling rate ----
            fs = abs(1/(snirf.data.time(2)-snirf.data.time(1)));
            
            % ---- Find trigger-based samples ----
            start_idx = find(tri_codes == params.trim.start_trigger, 1, 'first');
            %disp(start_idx);
            end_idx   = find(tri_codes == params.trim.end_trigger,   1, 'last');
            %disp(end_idx);

            if isempty(start_idx) || isempty(end_idx)
                error('Start or end trigger not found');
            end
            
            buffer_samples = round(params.trim.buffer*fs); 
            %disp(buffer)
            firstSample = tri_samples(start_idx) - buffer_samples;
            lastSample  = tri_samples(end_idx);
                   
            fprintf('Trim start: %d\n', firstSample);
            fprintf('Trim end:   %d\n', lastSample);

            % ---- Trim main data ----
            %isequal(snirf.data, snirf_trimmed.data)

            %fprintf('Size of original snirf before trimming data %d \n', ...
                %size(snirf.data.dataTimeSeries));

            snirf_trimmed.data.dataTimeSeries = ...
                snirf.data.dataTimeSeries(firstSample:lastSample, :);

            %isequal(snirf.data, snirf_trimmed.data)

            %fprintf('Size of original snirf after trimming data %d \n', ...
                %size(snirf.data.dataTimeSeries));

            %fprintf('Size of trimmed snirf after trimming data %d \n', ...
                %size(snirf_trimmed.data.dataTimeSeries));

            snirf_trimmed.data.time = ...
                snirf.data.time(firstSample:lastSample);

            %fprintf('Size of original snirf time after trimming %d \n', ...
                %size(snirf.data.time));

            %fprintf('Size of trimmed snirf time after trimming %d \n', ...
                %size(snirf_trimmed.data.time));

            % ---- Trim aux data ----
            if ~isempty(snirf.aux)
                for a = 1:length(snirf.aux)
                    if size(snirf.aux(a).dataTimeSeries, 1) >= lastSample
                        snirf_trimmed.aux(a).dataTimeSeries = ...
                            snirf.aux(a).dataTimeSeries(firstSample:lastSample, :);
                        %fprintf('Size of original aux after trimming data %d \n', ...
                            %size(snirf.aux(a).dataTimeSeries));
            
                        %fprintf('Size of trimmed aux after trimming data %d \n', ...
                            %size(snirf_trimmed.aux(a).dataTimeSeries));
                        
                        snirf_trimmed.aux(a).time = ...
                            snirf.aux(a).time(firstSample:lastSample, :);
                        %fprintf('Size of original aux time after trimming %d \n', ...
                            %size(snirf.aux(a).time));
            
                        %fprintf('Size of trimmed aux time after trimming %d \n', ...
                            %size(snirf_trimmed.aux(a).time));
                    end
                end
            end

            %% ===============================
            %% Raw Data Description
            %% ===============================

            % -------------------------------
            % 1) Number of channels
            % -------------------------------
            d = snirf_trimmed.data.dataTimeSeries;     % time × measurements
            t = snirf_trimmed.data.time; % time

            numchannels = size(d,2)/2; % number of channels

            % Channel labels
            chanLabels = arrayfun(@(i) ...
                sprintf('S%d-D%d', chanSrcDet(i,1), chanSrcDet(i,2)), ...
                1:numchannels, 'UniformOutput', false);

            % Get prope info 
            SD = struct();
            SD.MeasList = snirf_trimmed.data.cache.measurementListMatrix;
            SD.Lambda = snirf_trimmed.probe.wavelengths;
            SD.MeasListAct = SD.MeasList(:,3);
            SD.SrcPos = snirf_trimmed.probe.sourcePos2D;
            SD.DetPos = snirf_trimmed.probe.detectorPos2D;

            % -------------------------------
            % 2) Time elapsed & Sampling rate
            % -------------------------------
            duration_sec = t(end) - t(1); 
            duration_min = duration_sec / 60;
            fs = abs(1/(t(2)-t(1))); 
            
            % -------------------------------
            % 3) Short-Channel detection
            % -------------------------------
            idx_short = params.short_channels.map;

            % -------------------------------
            % 4) Accelerometer detection
            % -------------------------------
            numaux = numel(snirf_trimmed.aux);
            %fprintf('Aux present: %d entries\n', numaux);
            
            has_aux = numaux > 0;
            
            if has_aux
            
                % -------------------------------
                % Detect number of accelerometers
                % -------------------------------
                % Assume each accelerometer has 3 axes (x,y,z)
                n_axes = 3;
            
                n_accel = floor(numaux / n_axes);
                %fprintf('Detected %d accelerometer(s)\n', n_accel);
            
                t_fnirs = t(:);
                accelerometer_all = [];
            
                % -------------------------------
                % Loop over each accelerometer
                % -------------------------------
                for acc = 1:n_accel
            
                    % Channel indices for this accelerometer
                    ch_idx = (acc-1)*n_axes + (1:n_axes);
            
                    % Length from first channel
                    Nt = length(snirf_trimmed.aux(ch_idx(1)).dataTimeSeries);
            
                    % Preallocate
                    accel = zeros(Nt, n_axes);
            
                    % Load raw data
                    for ch = 1:n_axes
                        accel(:, ch) = snirf_trimmed.aux(ch_idx(ch)).dataTimeSeries;
                    end
            
                    % -------------------------------
                    % Interpolate to fNIRS time base
                    % -------------------------------
                    accel_i = zeros(length(t_fnirs), n_axes);
            
                    for ch = 1:n_axes
                        t_acc = snirf_trimmed.aux(ch_idx(ch)).time(:);
                        x_acc = snirf_trimmed.aux(ch_idx(ch)).dataTimeSeries(:);
            
                        accel_i(:, ch) = interp1( ...
                            t_acc, x_acc, t_fnirs, 'linear', 'extrap');
                    end
            
                    % -------------------------------
                    % Low-pass filter
                    % -------------------------------
                    FilterType = 1;
                    FilterOrder = 3;
                    [b,a] = MakeFilter(FilterType, FilterOrder, fs, 0.5, 'low');
            
                    accel_f = filtfilt(b, a, accel_i);
            
                    % -------------------------------
                    % Append to final matrix
                    % -------------------------------
                    accelerometer_all = [accelerometer_all, accel_f];
            
                end
            
                % Final output
                accelerometer = accelerometer_all;
            
            end

            % -------------------------------
            % 5) Print description summary
            % -------------------------------
            fprintf('\nRAW SUMMARY:\n');
            fprintf('----------------------------\n');
            fprintf('Number of optode channels: %d\n', numchannels);
            fprintf('Recording duration: %.2f minutes\n', duration_min);
            fprintf('Sampling rate: %.2f Hz\n', fs);
            
            if has_aux
                fprintf('Detected %d accelerometer(s)\n', n_accel);
            else
                fprintf('No accelerometer channels detected\n');
            end
            
            if ~isempty(idx_short)
                fprintf('Short-separation channels (optode index): %s\n', mat2str(idx_short));
            else
                fprintf('No short-separation channels detected\n');
            end
            fprintf('----------------------------\n');

            %% ===============================
            %% Channel Quality Assessment
            %% ===============================
            % Metrics:
            % 1) Detector Saturation
            % 2) Proportion of clean motion (motion artifact detection)
            % 3) SNR (dB)
            % 4) Scalp Coupling Index (SCI)
            % 5) Peak Spectral Power (PSP) 
            % 6) Quality Testing of Near-Infrared Scans <SCI & PS>

            % -------------------------------
            % Band-pass for cardiac signal
            % -------------------------------
            FilterType  = 1;   % Butterworth
            FilterOrder = 4;
            freq_range = get_thresh(params, "heart_band", [0.5, 2.5]);
            win_size = get_thresh(params, "window_size", 5);

            % Cardiac band: wide: 0.5-2.5 (Default); narrow: 0.7–1.5/2 Hz
            [card_b, card_a] = MakeFilter(FilterType, FilterOrder, fs, freq_range, 'bandpass');
            %[card_b, card_a] = MakeFilter(FilterType, FilterOrder, fs, 0.5, 'high');
            
            % Apply filters (zero-phase)
            dheart  = filtfilt(card_b, card_a, d);   % cardiac oscillations
            
            % -------------------------------
            % Preallocate metrics
            % -------------------------------
            quality_report = zeros(8,numchannels);
            SNR_linear      = zeros(numchannels,1);
            SNR_dB          = zeros(numchannels,1);
            saturated_detectors = false(numchannels,1);
            motion_clean    = zeros(numchannels,1); %motion artifacts
            % Select bad channel metric parameter
            selected_metrics = params.qc.metrics;
            % Select bad subject threshold parameter
            bad_sub_thresh  = get_thresh(params, "badsubject", 0.5);


            %% ===============================
            %% Converting Raw Light Intensity 
            %% to Optical Densities Changes
            % ===============================
            if params.nirs_use
                dod = hmrR_Intensity2OD_Nirs(d);
            else
                dod = hmrR_Intensity2OD(snirf_trimmed.data.copy);
                        end

            %% ===============================
            %% Motion Artifact Detection 
            %% ===============================
            tMotion = params.motion.detect.tMotion;   % detection window (seconds)
            tMask = params.motion.detect.tMask;     % temporal mask around artifact (seconds)
            STDEVthresh = params.motion.detect.STDEV;     % standard deviation threshold
            AMPthresh = params.motion.detect.AMP;       % amplitude threshold (optical density units)

           if params.nirs_use
                % Homer2-style
                tIncMan = ones(size(d,1),1);
                [tInc, tIncCh] = hmrMotionArtifactByChannel(d, fs, SD, tIncMan, ...
                    tMotion, tMask, STDEVthresh, AMPthresh);
            
                motion_map = tIncCh;
            
            else
                % Homer3-style
                [tIncAuto, tIncAutoCh] = hmrR_MotionArtifactByChannel( ...
                    dod, snirf_trimmed.probe, [], [], [], ...
                    tMotion, tMask, STDEVthresh, AMPthresh);
            
                motion_map = tIncAutoCh{1};
            end

            % -------------------------------
            % Compute metrics per channel
            % -------------------------------
            for c = 1:numchannels    

                % ---- 1) Detector saturation ----
                nanvec1 = isnan(d(:,c));
                nanvec2 = isnan(d(:,c + numchannels));
            
                if ~isempty(strfind(nanvec1', true(1, round(2*fs)))) || ...
                   ~isempty(strfind(nanvec2', true(1, round(2*fs))))
                   saturated_detectors(c) = true;
                else
                   saturated_detectors(c) = false;
                end
                quality_report(1,c) = saturated_detectors(c);

                % ---- 2) Proption of Clean Motion (Lack of Artifacts)  ----
                cp1 = mean(motion_map(:,c));                 % HbO clean %
                cp2 = mean(motion_map(:,c + numchannels));   % HbR clean %
                motion_clean(c) = mean([cp1 cp2]);
                quality_report(2,c) = motion_clean(c);

                % ---- 3) Signal to Noise Ratio (linear -> decibels [dB]) ----
                mu1 = mean(d(:,c));
                mu2 = mean(d(:,c + numchannels));
                sd1 = std(d(:,c));
                sd2 = std(d(:,c + numchannels));
            
                snr1 = mu1 / sd1;
                snr2 = mu2 / sd2;
            
                SNR_linear(c) = mean([snr1, snr2]);
                %SNR_dB(c) = 20 * log10(SNR_linear(c));

                snr1_dB = 20 * log10(snr1);
                snr2_dB = 20 * log10(snr2);
                
                SNR_dB(c) = mean([snr1_dB, snr2_dB]);
                
                quality_report(3,c) = SNR_dB(c);

                % QT-NIRS Cardiac QC:
                % 4) Scalp Coupling Index 
                % 5) Peak Spectral Power
                % Quality Testing score (proportion of SCI and PSP window
                % flags)

                % QT-NIRS-style thresholds
                sciThreshold = get_thresh(params, "SCI", 0.7);
                pspThreshold  = get_thresh(params, "PSP", 0.1);
                
                winSec  = win_size;
                winSamp = round(winSec * fs);
                nWin    = floor(size(dheart,1) / winSamp);
                
                QTscore  = nan(1, numchannels);

                x1 = dheart(:, c);
                x2 = dheart(:, c + numchannels);
                

                x1 = x1 ./ std(x1);
                x2 = x2 ./ std(x2);

                SCIwin = nan(nWin, 1);
                PSPwin = nan(nWin, 1);
                
                for w = 1:nWin
                    idx1 = (w-1)*winSamp + 1;
                    idx2 = w*winSamp;
                
                    w1 = x1(idx1:idx2);
                    w2 = x2(idx1:idx2);
                
                    similarity = xcorr(w1, w2, 'unbiased');
                    
                    if any(abs(similarity) > eps)
                        similarity = length(w1) * similarity ./ ...
                            sqrt(sum(abs(w1).^2) * sum(abs(w2).^2));
                        similarity(isnan(similarity)) = 0;
                
                        SCIwin(w) = similarity(length(w1));
                        %R = corrcoef(w1, w2);
                        %SCIwin(w) = R(1, 2);
                        
                        [pxx, f] = periodogram(similarity, hamming(length(similarity)), ...
                            length(similarity), fs, 'power');
                         bandMask = (f >= 0.5) & (f <= 2.5);
                
                        if any(bandMask)
                            PSPwin(w) = max(pxx(bandMask));
                        else
                            PSPwin(w) = NaN;
                        end
                     end
                 end
                
                goodWin = (SCIwin >= sciThreshold) & (PSPwin >= pspThreshold);
                
                QTscore(c) = mean(goodWin, 'omitnan');
                
                quality_report(4,c) = mean(SCIwin, 'omitnan');
                quality_report(5,c) = mean(PSPwin, 'omitnan');
                quality_report(6,c) = QTscore(c);


                % ---- 7) Bad Channel Criterion ----   
                prop_bad_selected = false;
               
                
                % ---- LOOP THROUGH USER METRICS ----
                for m = 1:length(selected_metrics)
                
                    metric = selected_metrics(m);
                
                    switch metric

                        case "SNR"
                            thresh = get_thresh(params, "SNR", 20);
                            value = quality_report(3,c);
                            bad = value < thresh;

                        case "SCI"
                            thresh = get_thresh(params, "SCI", 0.7);
                            value = quality_report(4,c);
                            bad = value < thresh;
               
                        case "PSP"
                            thresh = get_thresh(params, "PSP", 0.1);
                            value = quality_report(5,c);
                            bad = value < thresh;

                        case "QT"
                            thresh = get_thresh(params, "QT", 0.75);
                            value = quality_report(6,c);
                            bad = value < thresh;
                
                        otherwise
                            warning('Unknown QC metric: %s', metric);
                            continue
                    end
                
                
                    % ---- Determine if selectived metric passed threshold ----
                    prop_bad_selected = prop_bad_selected | bad;

                end
                
                is_bad_chan = ...
                    saturated_detectors(c) || ...
                    (motion_clean(c) < get_thresh(params, "motion_clean", 0.6)) || ...
                    prop_bad_selected;
                
                quality_report(7,c) = is_bad_chan;
            end

            % ---- 8) Bad Subject Criterion ----
            
            fprintf('\nCHANNEL QUALITY SUMMARY:\n');
            fprintf('----------------------------\n');

            % ---- ALWAYS USED ----
            sub_saturated = any(quality_report(1,:) == 1);
            sub_bad_motion = mean(quality_report(2,:) < get_thresh(params, "motion_clean", 0.6), 'omitnan');
            
            % ---- Gather bad channels from subject ----
            bad_channels = quality_report(7,:) == 1;
            
            % --- Proportion of bad channels
            prop_bad = mean(bad_channels);
            
            fprintf('Bad Channels: %d (%.1f%%)\n', ...
                sum(bad_channels), 100 * prop_bad);
            
            if any(bad_channels)
                fprintf('Channels: %s\n', strjoin(chanLabels(bad_channels), ', '));
            end
            
            % ---- Bad subject ----
            is_bad_subject = prop_bad >= bad_sub_thresh;
            
            if is_bad_subject
                fprintf('Subject flagged as BAD (%.2f >= %.2f)\n', prop_bad, bad_sub_thresh);
            else
                fprintf('Subject is GOOD (%.2f < %.2f)\n', prop_bad, bad_sub_thresh);
            end
            
            quality_report(8,:) = is_bad_subject;

            quality_report = array2table(quality_report, ...
                'VariableNames', chanLabels, ...
                'RowNames', cellstr(metrics));
            
            quality_report.Properties.DimensionNames = ["Metric", "Channel"];

            all_quality_reports{end +1} = quality_report;
            all_subject_ids{end +1} = subjectID;

            % ---- COUNTERS ----
            subject_total_count = subject_total_count + 1;
            if is_bad_subject
                subject_fail_count = subject_fail_count + 1;
            end
            
            % ---- PRINT FINAL SUMMARY ----
            fprintf('Subject QC → Saturated=%d | Motion=%.2f | SelectedFlag=%d | Final=%d\n', ...
                sub_saturated, sub_bad_motion, prop_bad_selected, is_bad_subject);
                        fprintf('----------------------------\n');

            fprintf('\nPREPROCESS SUMMARY:\n');
            fprintf('----------------------------\n');
            %% ===============================
            %% Motion Artifact Correction 
            %% ===============================
            
            method = params.motion.correct.method;
            
            switch method
                
                % ===============================
                % Temporal Derviative Distribution Repair
                % ===============================
                case "tddr"
            
                    dod_corrected = TDDR(dod.dataTimeSeries, fs);
                    dodMC = dod;
                    dodMC.dataTimeSeries = dod_corrected;
                    fprintf('Motion Correction: Temporal Derviative Distribution Repair \n');

                % ===============================
                % Targetted PCA (Homer3)
                % ===============================
                case "tPCA"
            
                    nSV = get_motion_param(params, "tPCA", "nSV", 0.97);
                    maxIter = get_motion_param(params, "tPCA", "maxIter", 3);
           
                    dodMC = hmrR_MotionCorrectPCArecurse(dod.copy, snirf_trimmed.probe, {}, {}, {}, tMotion, tMask, STDEVthresh, AMPthresh, nSV, maxIter, 1);
%
                    fprintf('Motion Correction: Targetted PCA (nSV=%.2f)\n', nSV);

                % ===============================
                % rLOESS (Homer3)
                % ===============================
                case "rLOESS"
            
                    span = get_motion_param(params, "rLOESS", "span", 0.02);
           
                    dodMC = hmrR_MotionCorrectRLOESS(dod.copy, span, 1);
            
                    fprintf('Motion Correction: rLOESS (span=%.2f)\n', span);

                % ===============================
                % Wavelet (Homer3)
                % ===============================
                case "wavelet"
            
                    iqr = get_motion_param(params, "wavelet", "iqr", 1.5);
            
                    dodMC = hmrR_MotionCorrectWavelet(dod.copy, [], [], iqr, 1);
            
                    fprintf('Motion Correction: Wavelet (iqr=%.2f)\n', iqr);
            
            
                % ===============================
                % Spline + SG (Homer3)
                % ===============================
                case "splineSG"
            
                    p = get_motion_param(params, "splineSG", "p", 0.99);
                    FrameSize_sec = get_motion_param(params, "splineSG", "FrameSize_sec", 10);
            
                    dodMC = hmrR_MotionCorrectSplineSG(dod, {}, p, FrameSize_sec);
            
                    fprintf('Motion Correction: Homer3 SplineSG (p=%.2f, Frame=%.2f)\n', p, FrameSize_sec);
            
            
                % ===============================
                % Custom Spline + SG (Homer2-style)
                % ===============================
                case "splineSG_custom"
            
                    p = get_motion_param(params, "splineSG_custom", "p", 0.99);
                    FrameSize_sec = get_motion_param(params, "splineSG_custom", "FrameSize_sec", 2);
                    K = get_motion_param(params, "splineSG_custom", "K", 3);
          
                    % ---- spline correction ----
                    dspline = hmrR_MotionCorrectSpline_Nirs(dod, t, SD, tIncCh, p);
            
                    % ---- SG filter ----
                    FrameSize = round(FrameSize_sec * fs);
                    if mod(FrameSize,2)==0
                        FrameSize = FrameSize + 1;
                    end
            
                    dsg = sgolayfilt(dspline, K, FrameSize);
            
                    dodMC = dsg;
            
                    fprintf('Motion Correction: Custom Spline+SG (p=%.2f, Frame=%d, K=%d)\n', p, FrameSize, K);
            
            
                % ===============================
                % No correction
                % ===============================
                case "none"
            
                    dodMC = dod;
            
                    fprintf('Motion Correction: NONE\n');
            
            
                otherwise
                    error('Unknown motion correction method: %s', method);
            end
            
            if  params.motion.plot.enable
            
                % Raw OD
                if params.nirs_use
                    dRaw = dod;   % Homer2 style (matrix)
                else
                    dRaw = dod.dataTimeSeries;   % Homer3
                end
            
                % Motion-corrected OD
                if params.nirs_use
                    dMC = dodMC;
                else
                    dMC = dodMC.dataTimeSeries;
                end
            
                % Time
                t_plot = t;
            
                % Channels to plot
                if isfield(params.motion.plot, 'channels') && ~isempty(params.motion.plot.channels)
                    chIdx = params.motion.plot.channels;
                else
                    chIdx = 1:min(5, numchannels); % default
                end
            
                % Plot
                plotODChannels(t_plot, dRaw, dMC, chIdx, numchannels, motion_clean, ...
                    sprintf('%s — ΔOD Motion Correction', subjectID));
            end

            %% ===============================
            %% Optical Density Changes
            %% to Hemoglobin Concentration Changes
            %% ===============================

            % Computing the concentration changes in O2Hb/HbO (oxygenated hemoglobin) and HHb/HbR (deoxygenated hemoglobin)
            % by means of the Modified Lambert-Beer Law <https://iopscience.iop.org/article/10.1088/0031-9155/33/12/008 
            % (Delpy et al., 1988)>.  This conversion requires the user to insert ppf (partial 
            % pathlength factor) values corresponding to each wavelength. These are age and 
            % wavelength specific. Here we use 6.06.
        

            % Partial pathlength factor    
            ppf = params.hb.ppf;
            
            if strcmp(params.motion.correct.method, 'splineSG_custom')
                % Homer2-style
                dhb = hmrR_OD2Conc_Nirs(dodMC, SD, ppf);
            else
                % Homer3-style
                dhb = hmrR_OD2Conc(dodMC, snirf_trimmed.probe, ppf);
            end
            
            fprintf('Hb Conversion: PPF = [%g %g]\n', ppf(1), ppf(2));
            
            % ===============================
            % Bandpass Filtering 
            % ===============================
            
            hpf  = params.filter.hpf;
            lpf = params.filter.lpf;
            
            if params.filter.use
            
                if strcmp(params.motion.correct.method, 'splineSG_custom')
                    % Homer2 filtering
                    dhbFilt = hmrR_BandpassFilt_Nirs(dhb, fs, hpf, lpf);
                else
                    % Homer3 filtering
                    dhbFilt = hmrR_BandpassFilt(dhb, hpf, lpf);
                end
            
                fprintf('Filtering: Bandpass [%.3f %.3f]\n', hpf, lpf);
            else
                    % No filtering
                    dhbFilt = dhb;
                    fprintf('No Filtering');
            end
            fprintf('----------------------------\n');
            
            %% ===============================
            %% Extract Hemoglobin Types
            %% ===============================
            if strcmp(params.motion.correct.method, 'splineSG_custom')
                HbO = squeeze(dhbFilt(:,1,:));
                HbR = squeeze(dhbFilt(:,2,:));
                HbT = squeeze(dhbFilt(:,3,:));
            else    
                % filtered → always dataTimeSeries
                HbO = dhbFilt.dataTimeSeries(:, 1:numchannels);
                HbR = dhbFilt.dataTimeSeries(:, numchannels+1 : 2*numchannels);
                HbT = dhbFilt.dataTimeSeries(:, 2*numchannels+1 : 3*numchannels);
            end

            % ===============================
            % OPTIONAL Hb PLOTTING
            % ===============================
            if  params.hb.plot.enable
            
                % Select channels
                if ~isempty(params.hb.plot.channels)
                    chIdx = params.hb.plot.channels;
                else
                    chIdx = 1:min(5, numchannels);
                end
            
                % Get bad channels (logical)
                is_bad_chan = logical(table2array(quality_report("BadChan", :)));
            
                % ---- Plot 1: Good vs Bad ----
                plotHbGoodBad(t, HbO, HbR, chIdx, is_bad_chan, subjectID);
            
                % ---- Plot 2: Spaghetti + Mean ----
                plotHbSpaghetti(t, HbO, HbR, 1:numchannels, subjectID);
            
            end



            outfile_tmp = fullfile(outpath, ...
                [subjectID '_raw_segmented.mat']);

            save(outfile_tmp, ...
                    'snirf_trimmed');

            outfile = fullfile(outpath, ...
                [subjectID '_preprocessed.mat']);
            
            if exist(outfile, 'file')
                delete(outfile);
            end

            if has_aux
                save(outfile, ...
                    'HbO', 'HbR', 'HbT', 'dhbFilt','tri_raw', 'idx_short', ...
                    'accelerometer', 'quality_report', 'fs', 't', 'SD', 'probeInfo', '-v7.3');
            else
                save(outfile, ...
                    'HbO', 'HbR', 'HbT', 'dhbFilt', 'tri_raw', 'idx_short', ...
                    'quality_report', 'fs', 't', 'SD', 'probeInfo', '-v7.3');
            end
            
            end
        end
    end

%% ===============================
%% FINAL GROUP QC ASSESSMENT
%% ===============================
if ~isempty(all_quality_reports)
    quality_report_assessment(all_quality_reports, all_subject_ids, outdir, metrics);
else
    warning('No quality reports collected — skipping group QC summary');
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

function plotODChannels(t, dRaw, dMC, chIdx, numchannels, motion_clean, figTitle)

if isempty(chIdx)
    warning('%s: no channels to plot', figTitle);
    return
end

for c = chIdx

    % wavelength split
    wl1 = c;                  % 760 nm
    wl2 = c + numchannels;    % 850 nm

    % scale
    pad = 1.1;
    ylims = pad * max(abs([ ...
        dRaw(:,[wl1 wl2]); ...
        dMC(:,[wl1 wl2]) ]), [], 'all');

    figure('name', sprintf('%s | Ch %d', figTitle, c), 'color','w');

    % ================= RAW =================
    ax1 = subplot(2,1,1); hold on
    plot(t, dRaw(:,wl1), 'r')
    plot(t, dRaw(:,wl2), 'b')
    title(sprintf('RAW | Ch %d | MotionClean=%.2f', c, motion_clean(c)))
    ylabel('\DeltaOD')
    ylim([-ylims ylims])
    legend('760 nm','850 nm')

    % ================= MC =================
    ax2 = subplot(2,1,2); hold on
    plot(t, dMC(:,wl1), 'r')
    plot(t, dMC(:,wl2), 'b')
    title('Motion Corrected')
    xlabel('Time (s)')
    ylabel('\DeltaOD')
    ylim([-ylims ylims])

    linkaxes([ax1 ax2],'x')
    pause(0.5)
end
end

function plotHbGoodBad(t, HbO, HbR, chIdx, is_bad_chan, subjectID)

figure('Name', [subjectID ' | Hb Good vs Bad'], 'Color','w');

for c = chIdx

    subplot(length(chIdx),1,find(chIdx==c)); hold on

    if is_bad_chan(c)
        plot(t, HbO(:,c), 'r', 'LineWidth', 1.5);
        plot(t, HbR(:,c), 'm', 'LineWidth', 1.5);
        title(sprintf('Ch %d (BAD)', c))
    else
        plot(t, HbO(:,c), 'b');
        plot(t, HbR(:,c), 'c');
        title(sprintf('Ch %d (GOOD)', c))
    end

    ylabel('\DeltaHb')
end

xlabel('Time (s)')

% global legned (dummy variables)
hold on

hHbO_good = plot(nan, nan, 'b');
hHbR_good = plot(nan, nan, 'c');
hHbO_bad  = plot(nan, nan, 'r', 'LineWidth', 1.5);
hHbR_bad  = plot(nan, nan, 'm', 'LineWidth', 1.5);

legend([hHbO_good, hHbR_good, hHbO_bad, hHbR_bad], ...
       {'HbO (good)', 'HbR (good)', 'HbO (bad)', 'HbR (bad)'}, ...
       'Position', [0.85 0.9 0.1 0.1])  % adjust as needed

end

function plotHbSpaghetti(t, HbO, HbR, chIdx, subjectID)

figure('Name', [subjectID ' | Hb Spaghetti'], 'Color','w');

% ---- HbO ----
subplot(2,1,1); hold on
for c = chIdx
    plot(t, HbO(:,c), 'Color', [0 0 1 0.1]) % transparent blue
end

meanHbO = mean(HbO(:,chIdx),2,'omitnan');
plot(t, meanHbO, 'b', 'LineWidth', 3)

title('HbO (Spaghetti + Mean)')
ylabel('\DeltaHb')

% ---- HbR ----
subplot(2,1,2); hold on
for c = chIdx
    plot(t, HbR(:,c), 'Color', [1 0 0 0.1]) % transparent red
end

meanHbR = mean(HbR(:,chIdx),2,'omitnan');
plot(t, meanHbR, 'r', 'LineWidth', 3)

title('HbR (Spaghetti + Mean)')
xlabel('Time (s)')
ylabel('\DeltaHb')

end

function quality_report_assessment(all_quality_reports, subject_ids, outdir, metrics)

nSubj = length(all_quality_reports);

fprintf('\n================ GROUP QC SUMMARY ================\n');


% --- Subject-level summary --- 
bad_subject_flags = zeros(nSubj,1);

for s = 1:nSubj
    qr = all_quality_reports{s};
    bad_subject_flags(s) = any(table2array(qr("BadSub", :)));
end

prop_bad_subjects = mean(bad_subject_flags);

fprintf('Bad Subjects: %d / %d (%.2f%%)\n', ...
    sum(bad_subject_flags), nSubj, 100*prop_bad_subjects);

% --- Metric-level summary --- 
fprintf('\nMETRIC SUMMARY:\n');

metric_names = [];
metric_values = [];

for m = 1:length(metrics)

    metric = metrics(m);

    if metric == "BadSub"
        continue
    end

    vals = [];

    for s = 1:nSubj
        qr = all_quality_reports{s};
        v = table2array(qr(metric, :));
        vals = [vals; v(:)];
    end

    if metric == "BadChan" || metric == "DetectorSat"
        val = mean(vals == 1);
        fprintf('%s bad proportion: %.2f%%\n', metric, 100*val);
    else
        val = mean(vals, 'omitnan');
        sd  = std(vals, 'omitnan');
        fprintf('%s mean: %.3f ± %.3f\n', metric, val, sd);
    end

    metric_names = [metric_names; metric];
    metric_values = [metric_values; val];
end

% ===============================
% CHANNEL-LEVEL SUMMARY
% ===============================
fprintf('\nCHANNEL-LEVEL SUMMARY:\n');

qr0 = all_quality_reports{1};
chanLabels = qr0.Properties.VariableNames;
nChannels = length(chanLabels);

channel_bad_matrix = zeros(nSubj, nChannels);

for s = 1:nSubj
    qr = all_quality_reports{s};
    channel_bad_matrix(s,:) = table2array(qr("BadChan", :));
end

channel_bad_rate = mean(channel_bad_matrix,1);

[sorted_vals, idx] = sort(channel_bad_rate, 'descend');

for i = 1:min(10, length(idx))
    fprintf('%s bad rate: %.2f%%\n', chanLabels{idx(i)}, 100*sorted_vals(i));
end

% --- Tabulating quality reports  --- 

% ---- 1. Channel-level table ----
T_channel = table(chanLabels', channel_bad_rate', ...
    'VariableNames', {'Channel', 'BadRate'});

% ---- 2. Metric-level table ----
T_metric = table(metric_names, metric_values, ...
    'VariableNames', {'Metric', 'Value'});

% ---- 3. Subject-level table ----
T_subject = table(subject_ids', bad_subject_flags, ...
    'VariableNames', {'SubjectID', 'IsBad'});

group_qc = struct();
group_qc.channel = T_channel;
group_qc.metric  = T_metric;
group_qc.subject = T_subject;

outfile_mat = fullfile(outdir, 'group_qc_summary.mat');
save(outfile_mat, 'group_qc');

fprintf(' - MATLAB QC object: %s\n', outfile_mat);
fprintf('=================================================\n');

end

% Please cite TDDR function if of use: 
%   Fishburn F.A., Ludlum R.S., Vaidya C.J., & Medvedev A.V. (2019). 
%   Temporal Derivative Distribution Repair (TDDR): A motion correction 
%   method for fNIRS. NeuroImage, 184, 171-179.
%   https://doi.org/10.1016/j.neuroimage.2018.09.025

function signal_corrected = TDDR( signal , sample_rate )
% This function is the reference implementation for the TDDR algorithm for 
%   motion correction of fNIRS data, as described in:
%
%   Fishburn F.A., Ludlum R.S., Vaidya C.J., & Medvedev A.V. (2019). 
%   Temporal Derivative Distribution Repair (TDDR): A motion correction 
%   method for fNIRS. NeuroImage, 184, 171-179.
%   https://doi.org/10.1016/j.neuroimage.2018.09.025
%
% Usage:
%   signals_corrected = TDDR( signals , sample_rate );
%
% Inputs:
%   signals: A [sample x channel] matrix of uncorrected optical density data
%   sample_rate: A scalar reflecting the rate of acquisition in Hz
%
% Outputs:
%   signals_corrected: A [sample x channel] matrix of corrected optical density data
%

%% Iterate over each channel
nch = size(signal,2);
if nch>1
    signal_corrected = zeros(size(signal));
    for ch = 1:nch
        signal_corrected(:,ch) = TDDR( signal(:,ch) , sample_rate );
    end
    return
end

%% Preprocess: Separate high and low frequencies
filter_cutoff = .5;
filter_order = 3;
Fc = filter_cutoff * 2/sample_rate;
signal_mean = mean(signal);
signal = signal - signal_mean;
if Fc<1
    [fb,fa] = butter(filter_order,Fc);
    signal_low = filtfilt(fb,fa,signal);
else
    signal_low = signal;
end
signal_high = signal - signal_low;

%% Initialize
tune = 4.685;
D = sqrt(eps(class(signal)));
mu = inf;
iter = 0;

%% Step 1. Compute temporal derivative of the signal
deriv = diff(signal_low);

%% Step 2. Initialize observation weights
w = ones(size(deriv));

%% Step 3. Iterative estimation of robust weights
while iter < 50
    
    iter = iter + 1;
    mu0 = mu;
    
    % Step 3a. Estimate weighted mean
    mu = sum( w .* deriv ) / sum( w );
    
    % Step 3b. Calculate absolute residuals of estimate
    dev = abs(deriv - mu);

    % Step 3c. Robust estimate of standard deviation of the residuals
    sigma = 1.4826 * median(dev);

    % Step 3d. Scale deviations by standard deviation and tuning parameter
    r = dev / (sigma * tune);
    
    % Step 3e. Calculate new weights accoring to Tukey's biweight function
    w = ((1 - r.^2) .* (r < 1)) .^ 2;

    % Step 3f. Terminate if new estimate is within machine-precision of old estimate
    if abs(mu-mu0) < D*max(abs(mu),abs(mu0))
        break;
    end

end

%% Step 4. Apply robust weights to centered derivative
new_deriv = w .* (deriv-mu);

%% Step 5. Integrate corrected derivative
signal_low_corrected = cumsum([0; new_deriv]);

%% Postprocess: Center the corrected signal
signal_low_corrected = signal_low_corrected - mean(signal_low_corrected);

%% Postprocess: Merge back with uncorrected high frequency component
signal_corrected = signal_low_corrected + signal_high + signal_mean;

end

function peakHz = plot_cardiac_peak_one_subject(d, t, bandLow, bandHigh)
% plot_cardiac_peak_one_subject
% d         : time x channels matrix
% t         : time vector (same length as rows of d)
% bandLow   : lower frequency bound for search (e.g., 0.5)
% bandHigh  : upper frequency bound for search (e.g., 2.5)

    if nargin < 3 || isempty(bandLow)
        bandLow = 0.5;
    end
    if nargin < 4 || isempty(bandHigh)
        bandHigh = 2.5;
    end

    % Ensure column vector
    t = t(:);

    % Estimate sampling rate from time vector
    dt = median(diff(t));
    fs = 1 / dt;

    fprintf('Estimated sampling rate: %.3f Hz\n', fs);

    % Remove mean from each channel
    x = detrend(d, 'linear');

    % Compute power spectrum for each channel
    nCh = size(x, 2);
    allPxx = [];
    fvec = [];

    for ch = 1:nCh
        [Pxx, f] = pwelch(x(:, ch), [], [], [], fs);

        if isempty(fvec)
            fvec = f;
            allPxx = zeros(length(Pxx), nCh);
        end

        allPxx(:, ch) = Pxx;
    end

    % Average spectrum across channels
    meanPxx = mean(allPxx, 2, 'omitnan');

    % Search for peak in cardiac band
    cardiacMask = (fvec >= bandLow) & (fvec <= bandHigh);
    fCard = fvec(cardiacMask);
    pCard = meanPxx(cardiacMask);

    [peakPower, idxPeak] = max(pCard);
    peakHz = fCard(idxPeak);

    fprintf('Dominant cardiac peak: %.3f Hz (%.1f bpm)\n', peakHz, peakHz * 60);
    fprintf('Peak power in band: %.6g\n', peakPower);

    % Plot
    figure;
    plot(fvec, 10*log10(meanPxx), 'LineWidth', 1.5);
    hold on;
    xline(bandLow, '--', 'Band low');
    xline(bandHigh, '--', 'Band high');
    xline(peakHz, 'r--', sprintf('Peak %.2f Hz', peakHz));
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title('Mean Power Spectrum Across Channels');
    grid on;
end
