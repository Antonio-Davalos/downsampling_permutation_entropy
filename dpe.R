" Downsampling Permutation Entropy "

" by Antonio Dávalos

Version 1.0 - 20200918

This script is a collection if Multiscale Permutation Entropy (MPE) [1] and 
Downsampling Permutation Entropy techniques. The main
functions are

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

These functions were performed in [3], using the Bearing Fault dataset provided
by [4].

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
  [4] Bechhoefer, E. Fault Data Sets. 
  Available online: https://www.mfpt.org/fault-data-sets/
  
"


#  ------------------------ Auxiliary Functions -----------------------------

# -----------
# Lehmer Code
# -----------

lehmer <- function(motif) {
  dmot <- ncol(motif)
  vfact = factorial( rep((dmot-1):1, (dmot-1):1) )
  
  idx_ref <- rep(1:(dmot-1),(dmot-1):1)
  
  idx <- rep(2:dmot,(dmot-1))
  idx <- idx[ as.vector( !upper.tri( matrix( rep(TRUE, (dmot-1)^2 ), nrow = (dmot-1) ) ) ) ]
  
  rt <- motif[,idx] < motif[,idx_ref]
  
  return( rt %*% vfact )
}

# -----------
# Pattern Segmentation
# -----------

patterns <- function(x,n,niv,nsig,d) {
  ones <- matrix(rep(TRUE,d^2), ncol = d)
  idx <- rbind( !upper.tri(ones), 
                matrix( rep(TRUE,((n-2*d)*d) ), ncol = d ) ,
                !lower.tri(ones) )
  ids <- sort( rep( seq(0,(nsig-1)*n,by = n), niv*d) )
  
  C <- rep( 1:n, d)
  C <- rep( C[as.vector(idx)], nsig ) + ids
  
  A <- array( x[C], c(niv,d,nsig))
  A <- matrix( aperm( A , c(1,3,2) ), nrow = niv*nsig, ncol = d)
  return(A)
}

# -----------
# Pattern Counts
# -----------

pcount <- function(B,nsig) {
  df <- max(B) + 1 
  ntot <- dim(B)[1]
  ct <- matrix(B, ncol = nsig )
  
  ynames <- 0:(df-1)
  ymat <- matrix( rep(ynames,(ntot/nsig)) , ncol = df, byrow=T )
  
  y <- lapply( 1:nsig, function(x) colSums( ymat == ct[,x] ) ) 
  y <- matrix( unlist(y, use.names = FALSE), nrow = nsig, byrow = TRUE )
  colnames(y) <- ynames
  
  return(y)
}

# -----------
# Moving Average Filter
# -----------
ma <- function(x, m){
  xma <- as.numeric(stats::filter(x, rep(1 / m, m), sides = 2))
  xma <- xma[!is.na(xma)]
  return(xma)
}

# -----------
# Downsampling
# -----------
dwnsamp <- function(x,tau){
  nivdwn <- floor(length(x)/tau)
  xdwn <- matrix( x[1:(nivdwn*tau)], ncol = tau, byrow = TRUE )
  return(xdwn)
}

# -----------
# Composite Coarse-graining
# -----------
ccg <- function(x, m){
  xma <- ma(x,m)             
  xcoarse <- dwnsamp(xma,m)
  return(xcoarse)
}

# -----------
# Shannon's Entropy
# -----------
shannon <- function(y,niv,dfact) {
  p <- y / niv
  l <- log(p)
  l[is.infinite(l)] <- 0
  H <- ( rowSums(-l*p) )  /(log(dfact))  
  return(H)
}


#  ------------------------ Main MPE Functions -----------------------

# -----------
# Permutation Entropy
# -----------
pe <- function(xin, d) {
  
  # Parameter settings
  n <- dim(xin)[1]
  nsig <- dim(xin)[2]
  niv <- n-d+1
  dfact <- factorial(d)
  x <- array( xin, c(n,nsig,1) )
  
  # Segmentation
  A <- patterns(x, n, niv, nsig, d)
  
  # Encoding the sequence
  B <- lehmer(A)
  y <- pcount(B,nsig)
  
  # Permutation Entropy
  H <- shannon(y,niv,dfact)
  
  return(H)
}

# -----------
# Permutation Entropy - for refine composite algorithms
# -----------
rcpe <- function(xin, d) {
  
  # Parameter settings  
  n = dim(xin)[1]
  nsig = dim(xin)[2]
  niv <- n-d+1
  dfact = factorial(d)
  xarr <- array( xin, c(n,nsig,1) )
  
  # Segmentation
  A <- patterns(xarr, n, niv, nsig, d)
  
  # Encoding the sequence
  B <- lehmer(A)
  y <- pcount(B,1)
  
  # Permutation Entropy
  H <- shannon(y,sum(y),dfact)
  
  return(H)
}

# -----------
# Multiscale Permutation Entropy
# -----------
mpe <- function(x, d, m){
  xcoarse = ccg(x, m)
  H = pe( matrix(xcoarse[,1], ncol = 1),d)
  return(H)
}

#  ------------------------ Refinements -----------------------

# -----------
# Composite Multiscale Permutation Entropy
# -----------
cmpe <- function(x, d, m){
  xcoarse = ccg(x, m)
  H <- pe(xcoarse,d)
  return(mean(H))
}

# -----------
# Refined Composite Multiscale Permutation Entropy
# -----------
rcmpe <- function(x, d, m){
  xcoarse = ccg(x, m)
  H <- rcpe(xcoarse, d)
  return(H)
}

#  ------------------------ Downsampling Entropy -----------------------

# -----------
# Composite Downsampling Permutation Entropy
# -----------
cdpe <- function(x, d, m){
  xdown = dwnsamp(x, m)
  H <- pe(xdown,d)
  return(mean(H))
}

# -----------
# Refined Composite Downsampling Permutation Entropy
# -----------
rcdpe <- function(x, d, m){
  xdown = dwnsamp(x, m)
  H <- rcpe(xdown, d)
  return(H)
}

#  ------------------------ Vectorized D -----------------------

vpe <- Vectorize(pe, vectorize.args = 'd', SIMPLIFY = TRUE,
                USE.NAMES = TRUE)
vrcpe <- Vectorize(rcpe, vectorize.args = 'd', SIMPLIFY = TRUE,
                USE.NAMES = TRUE)
