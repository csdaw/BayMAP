y.update.ran <-
function(alpha, nzone, z, TC){
  rtnorm(n=sum(nzone), mean=rep(alpha, nzone), sd = 1, lower=0)^(z[TC==1]==3) * rtnorm(n=sum(nzone), mean=rep(alpha, nzone), sd = 1, upper=0)^(1-(z[TC==1]==3))
}
