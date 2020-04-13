# MetabolicDifferentiation
Matlab 2017b supplemental files

<b>phaseBright.m</b>

Main file for phase bright quantification of mature spores. The script takes a tif (red,green,blue,phase) and returns results as an array composed of number of cells that completed sporulation, number of cells that are engulfing, total cells sporulating, as well as corresponding percentages. Figures displaying each thresholding step are intended for use in quality control.


<b>Engulfmdl.mat</b>

A classification model built in Matlab 2017b using a training set of engulfing spores. This model is used in phaseBright.mat as a filter to eliminate objects that are not actual spores.


<b>mc_fs_fluorescence_ratio.m</b>

Main file to measure GFP intensity ratio between forespores and mothercells. The script takes a tif file (red,green,blue,phase) and returns are array of measured mothercell intensity, forespore intensity and calculated ratio for each forespore/mothercell pair. Figures can be generated displaying each individual pair or all pairs numbered.


<b>filterSnakes1.mat/filterSnakes2.mat</b>

Two different classification models used to filter out false positive forespore/mothercell pairs in mc_fs_fluorescence_ratio.mat


<b>motherFit.mat</b>

Classification model for the identification of mothercell objects post image segmentation.


<b>sporeFit.mat</b>

Classification model for the identification of forespore objects post image segmentation.

<b>numberToGene</b>

Related NCBI accession number to gene name and groups them according to gene lists in \regulon_libraries\. Example formats for inputs found in example-protein_result.xlsx and example-sequence.xlsx.

