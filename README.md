## fNIRS-preprocess-analysis
USC SoCo lab preprocessing function for fNIRS data, collected with NIRx

# About
This repository provides a flexible, end-to-end fNIRS preprocessing and postprocessing pipeline built in MATLAB. It is designed for SNIRF-based datasets (e.g., NIRx) and integrates standard methods from Homer2/Homer3 along with custom quality control, motion correction, and regression workflows. Includes glm modelling as well. 

The pipeline is parameter-driven, meaning all behavior is controlled through a single configuration file:
👉 create_params()
