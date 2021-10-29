# dpe
Downsampling Permutation Entropy

by Antonio Dávalos

Version 1.0 - 20200918

This script is a collection if Multiscale Permutation Entropy (MPE)[1], and Downsampling Permutation Entropy techniques. The following functions are implemented in R,

  - MPE: Original Multiscale Permutation Entropy [1]
  - cMPE: Composite MPE [2]
  - rcMPE: Refined Composite MPE [2]
  - cDPE: Composite Downsampling Permutation Entropy [3]
  - rcDPE: Refined Composite DPE [3]. 

Inputs:

  - x: (numeric) array, which represent the sequential data of a signal or time
  series.
  - d: (positive integer) embedded dimension.
  - m: (positive integer) time scale.
  
Outputs:

  - H: Permutation Entropy

These functions were performed in [3], using the Bearing Fault dataset provided by [4].

References:

  [1] Bandt C., Pompe B., 2002. Permutation entropy: a natural complexity measure for time series. Physical review letters, APS.
  
  [2] A. Humeau-Heurtier, C. Wu and S. Wu, Refined Composite Multiscale Permutation Entropy to Overcome Multiscale Permutation Entropy Length Dependence, in IEEE Signal Processing Letters, vol. 22, no. 12, pp. 2364-2367, Dec. 2015, doi: 10.1109/LSP.2015.2482603
  
  [3] Dávalos A, Jabloun M, Ravier P, Buttelli O. Improvement of Statistical Performance of Ordinal Multiscale Entropy Techniques Using Refined Composite Downsampling Permutation Entropy. Entropy. 2021; 23(1):30. https://doi.org/10.3390/e23010030
  
  [4] Bechhoefer, E. Fault Data Sets.  Available online: https://www.mfpt.org/fault-data-sets/
