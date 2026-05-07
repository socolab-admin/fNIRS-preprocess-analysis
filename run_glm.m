function run_glm(data_postprocessed)

nSubj = length(data_postprocessed);
%disp(nSubj)
for s = 1:nSubj

    data = data_postprocessed{s};
    %disp(data)
    ID = data.subjectID;
    %disp(ID)

    fprintf('Processing subject %s (%d/%d)\n', ID, s, nSubj);

    % ===============================
    % Bad channel masking for glm_hmrR function
    % ===============================
    is_bad_chan = table2array(data.quality_report("BadChan", :)) == 1;
    mlActAuto{1} = ~is_bad_chan(:); % contains bad channels to ignore
    %disp(mlActAuto)

    %% *GLM analysis*
    % After performing basic data preprocessing steps, we will proceed with the 
    % GLM analysis. GLM stands for General Linear Model. In GLM analysis, a reference 
    % function, which consists of a linear combination of regressors, is fitted to 
    % the time course of each channel. The design matrix for the GLM reference function 
    % is created by convolving the boxcar function (boxcar time course acquires values 
    % of 1 when the condition is _on _and 0 when the condition is_ off_) with the 
    % hemodynamic response function. 
    % 
    % We will use _hmrR_GLM _function to run the GLM. Below, we set the parameters 
    % for the function. Please refer to Homer3 function documentation to see an explixit 
    % explanation of these parameters 
            
    trange = [-2 30]; % tPre (pre-stimulus onset), tPost (duration of stimulus past stimulus onsent)
    glmSolveMethod = 1; % this specifies the GLM solution method to use
    % 1. use ordinary least squares (Ye et al (2009). NeuroImage, 44(2), 428?447.)
    % 2. use iterative weighted least squares (Barker, Aarabi, Huppert (2013). Biomedical optics express, 4(8), 1366?1379.)
     % Note that we suggest driftOrder=0 for this method as otherwise it can produce spurious results.
    idxBasis = 1; % this specifies the type of basis function to use for the HRF
    % 1. a consecutive sequence of gaussian functions
    % 2. a modified gamma function convolved with a square-wave of duration given by the stim marker.
       % The modified gamma function is (exp(1)*(t-tau).^2/sigma^2) .* exp(-(tHRF-tau).^2/sigma^2)
    % 3. a modified gamma function and its derivative convolved with a square-wave of duration given by the stim marker.
    % 4.  GAM function from 3dDeconvolve AFNI convolved with  a square-wave of duration given by the stim marker. (t/(p*q))^p * exp(p-t/q)
       % Defaults: p=8.6 q=0.547; The peak is at time p*q.  The FWHM is about 2.3*sqrt(p)*q.
    paramsBasis = [1.0, 1.0]; % Parameters for the basis function (chosen via idxBasis)
    % if idxBasis=1 [stdev step ~ ~ ~ ~] where stdev is the width of the gaussian and step is the temporal spacing between consecutive gaussians
    % if idxBasis=2. [tau sigma] applied to both HbO and HbR
    % or [tau1 sigma1 tau2 sigma2] recommended values [0.1 3.0 1.8 3.0] where the 1 (2) indicates the parameters for HbO (HbR).
    % if idxBasis=3 [tau sigma] applied to both HbO and HbR
    % or [tau1 sigma1 tau2 sigma2] recommended values [0.1 3.0 1.8 3.0] where the 1 (2) indicates the parameters for HbO (HbR).
    % if idxBasis=4 [p q] applied to both HbO and HbR
    % or [p1 q1 p2 q2] recommended values [0.1 3.0 1.8 3.0] where the 1 (2) indicates the parameters for HbO (HbR).

    %%
    % In addition, various other regressors can be added to the model. For example, 
    % Homer3 allows regressing Short Separation Channels (SSC), model drifts and auxiliary 
    % signals. To enable these regressors, we will specify some more parameters. Specifically, 
    % the drift order will be set to 3. Since the SSC were included in the set-up 
    % of an example measurement used here, we will set the SSC threshold to 1 cm distance.
            
    driftOrder = 3; % How many drift regressors to use (polynomial order)
    rhoSD_ssThresh = 1; % max distance for a short separation channel measurement. Set =0
    % if you do not want to regress the short separation measurements. Follows the static estimate procedure described in Gagnon et al (2011).
    % NeuroImage, 56(3), 1362?1371.
    flagNuisanceRMethod = 2; 
    % 0. if short separation regression is performed with the nearest short separation channel.
    % 1. if performed with the short separation channel with the greatest correlation.
    % 2. if performed with average of all short separation channels.
    % 3. uses tCCA regressors for nuisance regression, in Aaux, mapped by rcMap, provided by hmr_tCCA()
            
    % The specified parameters can now be used to run the GLM.
    
    dhbFilt = data.dhbFilt;

    stim = data.stim;
    probe = ProbeClass(data.SD);
    
    [data_yavg, data_yavgstd, nTrials, data_ynew, data_yresid, data_ysum2, beta_blks, yR_blks, hmrstats] = ...
        hmrR_GLM(dhbFilt, stim, probe, mlActAuto, [], [], [], trange, glmSolveMethod, idxBasis, paramsBasis,...
        rhoSD_ssThresh, flagNuisanceRMethod, driftOrder, []);
    
    %% PRINT Final Results
    fprintf('\n================ GLM RESULTS ================\n');

    %% -----------------------------
    % 1. Basic info
    %% -----------------------------
    disp('Beta labels (conditions):');
    disp(hmrstats.beta_label);
    
    disp('Number of trials per condition:');
    disp(nTrials);
    
    %% -----------------------------
    % 2. Extract stats
    %% -----------------------------
    tvals = hmrstats.tval;
    pvals = hmrstats.pval;
    ml = hmrstats.ml;
    
    % HbO = 1, HbR = 2
    tvals_HbO = squeeze(tvals(:,:,1));
    pvals_HbO = squeeze(pvals(:,:,1));
    
    tvals_HbR = squeeze(tvals(:,:,2));
    pvals_HbR = squeeze(pvals(:,:,2));
    
    %% -----------------------------
    % 3. Correlation fit (model vs data)
    %% -----------------------------
    disp('Correlation fit (R) between data and GLM model:');
    
    for iBlk = 1:length(yR_blks)
        fprintf('\n--- Block %d ---\n', iBlk);
        R = yR_blks{iBlk};  % (#channels x Hb)
    
        for ch = 1:size(R,1)
            fprintf('Ch %d | HbO R=%.3f | HbR R=%.3f\n', ...
                ch, R(ch,1), R(ch,2));
        end
    end
    
    %% -----------------------------
    % 4. Beta values (effect size)
    %% -----------------------------
    disp('Beta values (effect size):');
    
    for iBlk = 1:length(beta_blks)
        fprintf('\n--- Block %d ---\n', iBlk);
    
        beta = beta_blks{iBlk};  % (#coeff × Hb × channels × conditions)
    
        for ch = 1:size(beta,3)
            fprintf('Channel %d:\n', ch);
    
            for cond = 1:size(beta,4)
                fprintf('  Cond %d | HbO beta=%.4f | HbR beta=%.4f\n', ...
                    cond, beta(1,1,ch,cond), beta(1,2,ch,cond));
            end
        end
    end
    
    %% -----------------------------
    % 5. Channel-by-channel stats
    %% -----------------------------
    fprintf('\nChannel-wise statistics:\n');
    
    for ch = 1:size(pvals_HbO,2)
        fprintf(['Ch %d (Src %d - Det %d) | ' ...
                 'HbO: p=%.4f, t=%.2f | ' ...
                 'HbR: p=%.4f, t=%.2f\n'], ...
            ch, ml(ch,1), ml(ch,2), ...
            pvals_HbO(1,ch), tvals_HbO(1,ch), ...
            pvals_HbR(1,ch), tvals_HbR(1,ch));
    end
    
    %% -----------------------------
    % 6. Significant channels
    %% -----------------------------
    fprintf('\nSignificant channels (p < 0.05):\n');
    
    sig_idx = find(pvals_HbO(1,:) < 0.05);
    
    for i = sig_idx
        fprintf('Ch %d (Src %d - Det %d) | p=%.4f\n', ...
            i, ml(i,1), ml(i,2), pvals_HbO(1,i));
    end
    
    %% -----------------------------
    % 7. Residual summary
    %% -----------------------------
    disp('Residual variance (quick check):');
    
    for iBlk = 1:length(data_yresid)
        yres = data_yresid(iBlk).GetDataTimeSeries('reshape');
    
        var_res = var(yres(:));
        fprintf('Block %d residual variance: %.4f\n', iBlk, var_res);
    end
    
    fprintf('=============================================\n');

    %% Plot Final Results in a Graph
    % Displaying NIRS-data is typically done by plotting both the average changes 
    % of oxy- and deoxyhemoglobin in time. Here, we plot the resulting average * on 
    % HbO and HbR signals
    % Extract plot data from output of last function 
            
    % These are plot figure positioning parameters
    px      = .01;
    py      = .30;
    sx      = .30;
    sy      = .40;
    xstep   = .15;
            
    % Now plot the GLM processed block average 
    tAvg = data_yavg.time;
    yAvg = data_yavg.GetDataTimeSeries('reshape');
    nChannels = size(yAvg,3);
    channels = 1:nChannels;
    hbPlotNames = {'HbO ', 'HbR ', 'HbT '};
    hbType   = [1, 2, 3];
    for iHb = 1:length(hbType)
        figure('menubar','none', 'numbertitle','off', 'name',[hbPlotNames{iHb}, 'Concentration Block Average'], ...
            'units','normalized','position',[px, py, sx, sy]);
        px = px + xstep;
        for iCh = channels
            plot(tAvg, yAvg(:, hbType(iHb), channels(iCh)));
            hold on;
        end
        xlim([min(tAvg), max(tAvg)]);
        hold off;
        pause(1)
    end
    xlim([min(tAvg), max(tAvg)]);
    hold off;

end    
end
