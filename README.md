Fiber Photometry Calcium Signal Analysis Pipeline
This repository contains a set of MATLAB scripts for analyzing fiber photometry calcium imaging data. The pipeline covers bleach correction, event-related analysis (feeding and drinking), and window-based correlation analysis with body temperature and locomotion.

Overview
The scripts are designed to process CSV data containing calcium signals, reference signals, time stamps, and behavioral events. They provide:

Bleach correction using exponential or linear fitting.

Calculation of F/F₀ and ΔF/F.

Statistical characterization of calcium signals across multiple cells.

Classification of active vs. inactive cells.

Event-triggered averaging for feeding and drinking bouts.

Window-based correlation between calcium fluctuations and physiological/behavioral variables (temperature, locomotion).

File Descriptions
multi_cell_calcium_analysis.m
Processes multi-cell calcium data. Performs bleach correction (exponential full, exponential combined, or linear fallback), computes F/F₀, statistical metrics (mean, std, CV, kurtosis, skewness, dynamic range), classifies cells as active/inactive based on CV, and generates correlation matrices, heatmaps, network plots, and statistical visualizations.

feeding_event_analysis.m
Analyzes feeding events. Reads behavioral data, detects 0→1 transitions, merges consecutive events into bouts (gap threshold 60 s), performs per-event reference fitting, computes ΔF/F, and generates event-triggered averages, heatmaps, per-event plots, and full-session overview. Exports event data to CSV.

drinking_event_analysis.m
Similar to feeding analysis but for drinking (lick) events. Detects lick events from 500 Hz data, analyzes lick interval distribution, merges into bouts (gap threshold 15 s), performs per-event reference fitting, computes ΔF/F, and generates event-triggered averages, heatmaps, and full-session overview. Exports event data to CSV.

temperature_locomotion_correlation.m
Basic window-based correlation analysis. Reads body temperature and locomotion data (15 s bins), computes windowed statistics (calcium std, mean, CV, temperature slope/mean, locomotion total), and correlates them with calcium signals using Spearman correlation. Produces scatter plots and time-series plots.

temperature_locomotion_correlation_advanced.m
Advanced version of the above. Includes normalized coefficient of variation (CV_norm), window CV, and detailed correlation analysis between temperature/locomotion and multiple calcium metrics. Generates additional plots for CV_norm, window CV, and their relationships with temperature.

Requirements
MATLAB (R2020a or later recommended)

Statistics and Machine Learning Toolbox (corr, kurtosis, skewness, prctile, ksdensity)

Optimization Toolbox (lsqcurvefit)

Signal Processing Toolbox (movvar, fillmissing, interp1)

Image Processing Toolbox (optional, for imagesc – available in base MATLAB)

Curve Fitting Toolbox (optional, polyfit and polyval are in base MATLAB)

Usage
Set file paths and parameters
At the beginning of each script, modify:

filename: path to the calcium data CSV.

channel: channel index (0–3) to select signal and reference columns.

offset_sec, duration_sec: time range to analyze (default 0–7200 s).

For event analysis: eating_filename or drink_file, column index, sampling rate, etc.

For correlation analysis: combined_filename for temperature/locomotion data.

Output directories for CSV export.

Run the script
Execute the script in MATLAB. Figures will be generated and data may be exported to CSV.

Data format

Calcium data CSV: columns for signal, reference, and time (last column). First row is a header.

Feeding/drinking data: CSV or Excel with a column of 0/1 events.

Temperature/locomotion data: CSV with two columns (temperature, locomotion) and no header. Time is assumed to be in 15 s bins.

Outputs
Corrected F/F₀ traces and bleach curves.

Statistical tables (mean, std, CV, kurtosis, skewness, dynamic range).

Activity classification (active/inactive cells based on CV threshold).

Correlation matrices, clustered heatmaps, and network graphs for multi-cell data.

Event-triggered averages, heatmaps, and per-event traces for feeding/drinking.

Window-based correlation scatter plots and time-series plots.

Exported CSV files containing event-aligned data (dFF_percent, behavior, etc.).

Notes
All scripts assume time in seconds and sampling rates as specified in the code.

Reference channel outlier repair is performed by thresholding and linear interpolation.

For multi-cell data, each cell is processed independently.

Ensure all file paths are correctly set before running.

The scripts use caxis; if you are using MATLAB R2022a or later, consider replacing with clim.

License
This project is for research purposes. Please cite appropriately if used in publications.

Contact
For questions or issues, please contact the author.
