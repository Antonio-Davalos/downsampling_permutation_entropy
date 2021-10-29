" Downsampling Permutation Entropy - Example "

" by Antonio Dávalos

Version 1.0 - 20200630

This script is a simple example to test and use the functions

  - mpe: Original Multiscale Permutation Entropy [1]
  - cmpe: Composite MPE [2]
  - rcmpe: Refined Composite MPE [2]
  - cdpe: Composite Downsampling Permutation Entropy [3]
  - rcdpe: Refined Composite DPE [3]. 

contained in the file dpe.R.

The script generates an array of size N (random gaussian noise or an arima time
series), and proceeds to obtain the multiscale permutation entorpy for the five
methods above.

References:
  [1] Bandt C., Pompe B., 2002. Permutation entropy: a natural complexity 
  measure for time series. Physical review letters, APS
  [2] A. Humeau-Heurtier, C. Wu and S. Wu, Refined Composite Multiscale
  Permutation Entropy to Overcome Multiscale Permutation Entropy Length
  Dependence, in IEEE Signal Processing Letters, vol. 22, no. 12, pp. 2364-2367,
  Dec. 2015, doi: 10.1109/LSP.2015.2482603.
  [3] Dávalos A, Jabloun M, Ravier P, Buttelli O. Improvement of Statistical
  Performance of Ordinal Multiscale Entropy Techniques Using Refined Composite
  Downsampling Permutation Entropy. Entropy. 2021; 23(1):30.
  https://doi.org/10.3390/e23010030

"

#---------------------------------
# Sources
source("dpe.R")

#---------------------------------
# Random vector generator

N = 100000
#x <- rnorm(N , mean=0, sd=1)     # Variance is not important here
x <- as.numeric(arima.sim(list(order = c(1,0,0), ar = 0.9), n = N))

#----------------------------------
# Entropy parameter setting

mmax=30   # Time Scale
d=5   # Embedding Dimension

#----------------------------------
# Multiscale Permutation Entropy 

#H0 = pe(x,d)
H = numeric(mmax)
Hc = numeric(mmax)
Hrc = numeric(mmax)
Hcd = numeric(mmax)
Hrcd = numeric(mmax)

# MPE
pmtot <- proc.time()
for (ii in 1:mmax){
  H[ii] = mpe(x,d,ii)
}
pmtot <- proc.time() - pmtot
print(pmtot["elapsed"])

# cMPE
pmtot <- proc.time()
for (ii in 1:mmax){
  Hc[ii] = cmpe(x,d,ii)
}
pmtot <- proc.time() - pmtot
print(pmtot["elapsed"])

# cDPE
pmtot <- proc.time()
for (ii in 1:mmax){
  Hcd[ii] = cdpe(x,d,ii)
}
pmtot <- proc.time() - pmtot
print(pmtot["elapsed"])

# rcMPE
pmtot <- proc.time()
for (ii in 1:mmax){
  Hrc[ii] = rcmpe(x,d,ii)
}
pmtot <- proc.time() - pmtot
print(pmtot["elapsed"])

# rcDPE
pmtot <- proc.time()
for (ii in 1:mmax){
  Hrcd[ii] = rcdpe(x,d,ii)
}
pmtot <- proc.time() - pmtot
print(pmtot["elapsed"])
