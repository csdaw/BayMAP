alpha.update <-
function(y, zone, nzone, J, beta, X, tau){
  ysum <- tapply(y, factor(zone), sum)[rank(unique(zone))]
  alpha.hat <- (c(X%*%beta)/tau^2 + ysum)/(1/tau^2 + nzone)
  V.theta <- 1/(1/tau^2 + nzone)
  rnorm(J, alpha.hat, sqrt(V.theta))
}
