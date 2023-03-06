
# Process Package README

## main image processing scripts.

### Description: 
* process raw images and output individual spot locations and barcoding info for DNA seqFISH+ and RNA seqFISH+ and also output sequential immunofluorescence intensity information at each DNA seqFISH+ spot. DNA seqFISH+ spots were further decoded by the script under "paint_decoding" folder.

Main scripts used for the processing are named as "example_Script_seqfish_pipeline_".

### Dependencies
1. Matlab Version R2019b
2. Fiji installed
	* can be downloaded at https://imagej.net/Fiji/Downloads
3. Segmentation in 2d (RoiSet.zip) or 3d (labeled images) generated from Cellpose (https://www.nature.com/articles/s41592-020-01018-x).
4. radialcenter.m by Dr. Parthasarathy (https://pages.uoregon.edu/raghu/particle_tracking.html).
