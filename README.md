# fNIRS-preprocess-analysis
USC SoCo lab preprocessing function for fNIRS data, collected with NIRx

## About
This repository provides a flexible, end-to-end fNIRS preprocessing and postprocessing pipeline built in MATLAB. It is designed for SNIRF-based datasets (e.g., NIRx) and integrates standard methods from Homer2/Homer3 along with custom quality control, motion correction, and regression workflows. Includes glm modelling as well. 

The pipeline is parameter-driven, meaning all behavior is controlled through a single configuration file:
👉 create_params()

## Start Here
Before running anything:

👉 1. Read create_params() carefully
👉 2. Adjust parameters for your dataset and study design

This file defines:
- Data paths and folder structure
- Quality control thresholds
- Motion detection & correction methods
- Filtering and hemoglobin conversion
- Postprocessing (nuissance regression, epoching, etc.)

## Pipeline Overview

* Preprocessing (run_preprocessing)
- Load SNIRF + trigger (.tri) files
- Trim recordings using triggers
- Compute channel quality metrics:
  SNR, SCI, motion, saturation, etc.
- Convert → OD
- Detect motion artifacts
- Apply motion correction (tPCA, wavelet, splineSG, etc.)
- Convert → HbO / HbR / HbT
- Apply bandpass filtering
- Save subject-level outputs
- Generate group-level QC summary

* Postprocessing (run_postprocessing)
- Skip bad subjects (optional)
- Mask bad channels
- Build stimulus design matrix
- Nuisance regression (ssc, accelerometers):
- Z-score normalization 
- Epoch extraction 
- Save outputs + per-epoch CSV files

## Data Structure
Raw data structure:
rawdir/
  G01/
    S01_*/
      *.snirf
      *.tri
      probeInfo.mat

Output structure:
outdir/
  G01/
    S01_*/
      S01_preprocessed.mat
      S01_postprocessed.mat
      epoch CSV files
  group_qc_summary.mat
