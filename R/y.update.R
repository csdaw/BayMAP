y.update <-
function(X, beta, z, TC){
  rtnorm(n=dim(X)[1], mean=X%*%beta, sd = 1, lower=0)^(z[TC==1]==3) * rtnorm(n=dim(X)[1], mean=X%*%beta, sd = 1, upper=0)^(1-(z[TC==1]==3))
}
