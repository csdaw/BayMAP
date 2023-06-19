beta.update.ran <-
function(alpha, tau, X){
  v <- ginv(t(X)%*%X + (tau^2)/(10^6)*diag(ncol(X)))
  c(rmvnorm(1, v %*% t(X) %*% alpha, tau^2 * v))
}
