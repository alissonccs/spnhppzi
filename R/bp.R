# Computes the Bernstein polynomial's bases. (note: for computation stability, b is not divided by zeta here)
bp <- function(time, max_time_id, degree, zeta,N,n) {
  # n <- length(time)
  y <- time/zeta
  y1 <-max_time_id/zeta
  b <- matrix(nrow=N, ncol=degree)
  B <- matrix(nrow=n, ncol=degree)
  for(k in 1:degree)
  {
    # b[,k] <- stats::dbeta(y, k, degree - k + 1)/zeta
    b[,k] <- stats::dbeta(y, k, degree - k + 1)
    B[,k] <- stats::pbeta(y1, k, degree - k + 1)
  }
  return(list(b=b, B=B))
}
