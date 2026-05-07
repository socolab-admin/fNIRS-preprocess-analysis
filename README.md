# fNIRS-preprocess-analysis  
USC SoCo lab preprocessing function for fNIRS data (NIRx / SNIRF)

---

## About
This repository provides a flexible, end-to-end **fNIRS preprocessing and postprocessing pipeline** built in MATLAB. It is designed for **SNIRF-based datasets** (e.g., NIRx NIRSport2) and integrates standard methods from **Homer2/Homer3** along with custom quality control, motion correction, regression workflows, and **GLM modeling**.

The pipeline is **parameter-driven**, meaning all behavior is controlled through a single configuration file:

👉 **`create_params()`**

---

## 🚨 Start Here
Before running anything:

👉 **1. Read `create_params()` carefully**  
👉 **2. Adjust parameters for your dataset and study design**

This file defines:
- Data paths and folder structure  
- Quality control thresholds  
- Motion detection & correction methods  
- Filtering and hemoglobin conversion  
- Postprocessing (nuisance regression, epoching, etc.)

⚠️ The pipeline depends heavily on correct parameterization.

---

## Pipeline Overview

### Preprocessing (`run_preprocessing`)
- Load SNIRF + trigger (`.tri`) files  
- Trim recordings using triggers  
- Compute channel quality metrics:  
  - SNR, SCI, motion, saturation, etc.  
- Convert → Optical Density (OD)  
- Detect motion artifacts  
- Apply motion correction (tPCA, wavelet, splineSG, rLOESS, etc.)  
- Convert → HbO / HbR / HbT  
- Apply bandpass filtering  
- Save subject-level outputs  
- Generate group-level QC summary  

---

### Postprocessing (`run_postprocessing`)
- Skip bad subjects (optional)  
- Mask bad channels  
- Build stimulus design matrix  
- Nuisance regression:
  - Short-channel regression (SSC)  
  - Accelerometer regression (optional)  
- Z-score normalization  
- Epoch extraction (with HRF delay)  
- Save outputs + per-epoch CSV files  

---

## Data Structure

### 📥 Raw Data
rawdir/
G01/
S01_*/
*.snirf
*.tri
probeInfo.mat


- Works with **SNIRF files** (e.g., NIRx NIRSport2 systems)
- `.tri` = trigger file for trimming + epoching  
- `probeInfo.mat` = source-detector mapping  

---

### 📤 Output Data
outdir/
G01/
S01_/
S01_preprocessed.mat
S01_postprocessed.mat
S01_EpochXX*.csv
group_qc_summary.mat


- `*_preprocessed.mat` → cleaned + filtered signals  
- `*_postprocessed.mat` → regression + design matrix  
- `EpochXX.csv` → condition-based epochs  
- `group_qc_summary.mat` → dataset-level QC metrics  

---

