# Solar Active Region Analysis Pipeline

This repository contains a MATLAB script for processing solar active region (AR) data using EUV observations and NSO/GONG measurements. The script performs tasks by automatically identifying solar ARs, calculating their properties, and visualizing results.

Features
- Active Region Detection: Uses thresholding and edge-detection techniques to identify ARs from input data.
- EUV AR Properties Extraction: Computes latitude, longitude, area, brightness, and tilt angle for each detected AR.
- GONG AR Properties Extraction: Computes latitude, longitude, area, and helioseismic phaseshift properties for each detected AR.
- Data Visualization: Generates visualizations of processed data with overlays for active regions.

Output Storage: Saves processed data as .fits images, .csv tables, and .png images.

Requirements
- MATLAB: Version R2018b or later is recommended.
- Toolboxes: Ensure the Image Processing Toolbox is installed.

Input Data
- map_304.mat: Colormap data.
- TH.mat: Threshold data in the format [yyyymmdd tt TH].
- seg_coord.mat: Region segmentation coordinates [yyyymmdd lat1 lat2 long1 long2].
- STEREO .fts files organized in the fts_27 directory.
- GONG/fqm .fts files.
