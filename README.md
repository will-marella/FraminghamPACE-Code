# FraminghamPACE Analysis Code

This repository contains the code used to construct and evaluate **FraminghamPACE**, as described in:

> Marella WT, Ryan CP, Corcoran D, Eckstein Indik C, Furuya A, Kobor MS, Sugden K, Caspi A, Moffitt TE, and Belsky DW. *An epigenetic speedometer to measure Pace of Aging: FraminghamPACE.* medRxiv. 2026. https://doi.org/10.64898/2026.07.07.26357388

## Repository structure

The code is organized approximately in the order used in the study:

* `01_FHSBiomarkers/`: Framingham Heart Study biomarker processing
* `02_PaceOfAgingConstruction/`: construction of the longitudinal Pace of Aging measure
* `03_FraminghamPACEConstruction/`: training of the DNA-methylation algorithms
* `04_Analysis/`: analyses reported in the manuscript and supplement

## Calculating FraminghamPACE

This repository contains the study construction and analysis code. It is **not** the package for calculating FraminghamPACE in new DNA-methylation datasets.

For installation and usage instructions, see the FraminghamPACE R package:

https://github.com/will-marella/FraminghamPACE
