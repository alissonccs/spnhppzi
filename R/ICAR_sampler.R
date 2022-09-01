## Create ICAR samples
## W adajency matrix with diagonal entries 0
## variance of the ICAR
sample_icar <- function(W,sig=1){
  
  num <- rowSums(W)
  n <- ncol(W)
  
  Q <- -W
  diag(Q) <- num
  
  Q_aux=eigen(Q)$vectors[,order(eigen(Q)$values)]
  
  D_aux=sort(eigen(Q)$values)
  
  rnd <- rnorm(n-1,0,sqrt(sig*(1/D_aux[-1])))
  rnd <- Q_aux%*%c(0,rnd)
  return(as.vector(rnd))
}
