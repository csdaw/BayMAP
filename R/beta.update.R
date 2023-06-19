beta.update <-
function(beta, X, z, y){
  v <- ginv(t(X)%*%X + 1/(10^6)*diag(ncol(X)))
  beta <- mvrnorm(n=1,mu= v %*% t(X)%*%y, Sigma=v)
  l <- list(est = beta)
  
  return(l)
}
