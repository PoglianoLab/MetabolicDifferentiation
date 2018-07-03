# MetabolicDifferentiation
Matlab 2017b supplemental files

<b>PhaseBright.m</b>

Main file for phase bright quantification of mature spores. The script takes an image file name as input (baseFileName) and returns results as an array composed of number of cells that completed sporulation, number of cells that are engulfing, total cells sporulating, as well as corresponding percentages. Figures displaying each thresholding step are intended for use in quality control.


<b>Engulfmdl2.mat</b>

A classification model built in Matlab 2017b using a training set of engulfing spores. This model is used in PhaseBright.mat as a filter to eliminate objects that are not actual spores.


<b>GFPRatio.m</b>

Main file to measure GFP intensity ratio between forespores and mothercells. The script takes an image file name as input (baseFileName) and returns are array of measured mothercell intensity, forespore intensity and calculated ratio for each forespore/mothercell pair. Figures can be generated displaying each individual pair or all pairs numbered.


<b>filterSnakesV1.mat/filterSnakesV2.mat</b>

Two different classification models used to filter out false positive forespore/mothercell pairs in GFPRatio.mat


<b>motherFitV6.mat</b>

Classification model for the identification of mothercell objects post image segmentation.


<b>sporeFitV4.mat</b>

Classification model for the identification of forespore objects post image segmentation.
