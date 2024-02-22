# DEFRT EEG Study Preprocessing Pipeline

## Overview
This preprocessing pipeline is designed for the DEEFRT task within the DEFRT EEG Study. It converts raw EEG data from MFF format to Matlab usable data, processes the data through various preprocessing steps, and prepares it for further analysis. The pipeline includes parsing photodiode signals to identify task events, filtering, downsampling, artifact subspace reconstruction (ASR), channel interpolation, average re-referencing, and Independent Component Analysis (ICA) for artifact removal.

## Requirements
- MATLAB (tested with MATLAB R2023a)
- EEGLAB (version 2023.1 recommended)
- ASR plugin for EEGLAB
- ICLabel plugin for EEGLAB (for component classification)

## Setup
Before running the pipeline, ensure that the required software and plugins are installed and that the paths in the script are correctly set to your local environment. The script assumes your data is stored in a specific directory structure. Adjust the `ROOT`, `rawEEG_dir`, and other directory variables as needed to fit your data organization.

## Running the Pipeline
The script can be executed as a function or run line-by-line for individual datasets. To run the script:

1. Open the script in MATLAB.
2. Adjust the `ROOT` variable to point to your study folder.
3. Verify that the `rawEEG_dir` and `PROC_EEG_DIR` paths are correct.
4. Run the script. If running as a loop, note that it will pause during the ICA artifact selection step for manual review.

### Key Steps
1. **Loading and Preparing EEG Data**: Converts raw MFF files to EEGLAB's data format.
2. **Event Handling**: Uses photodiode signals to mark relevant task events.
3. **Filtering**: Applies a 1 Hz highpass and a 58 Hz lowpass filter, then downsamples to 200 Hz.
4. **Artifact Subspace Reconstruction (ASR)**: Removes bad channels and corrects bursty artifacts. Removed channels are then interpolated.
5. **Independent Component Analysis (ICA)**: Identifies and allows for the removal of artifactual components. Includes an automatic labeling step via ICLabel and a manual review step.

## Configuration
- The `runManualICA_Rejection` variable controls whether ICA component rejection is manual (1) or skipped (0).
- The pipeline uses default parameters for filtering and downsampling, which can be adjusted according to your study's needs.

## Notes
- The pipeline includes error checks for missing directories and files, ensuring smooth execution.
- The script is set up to handle multiple EEG datasets sequentially but can be easily adapted for single dataset processing.
- During manual ICA component rejection, the script will pause and wait for user input in the MATLAB Command Window.

## Updates
- **6/22/23**: Updated DEFRT script for alignment with GNG script.
- **6/27/23**: Corrected ASR step to disable window criterion.
