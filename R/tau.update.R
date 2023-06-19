tau.update <-
function(alpha, beta, X, J){
  sqrt(sum((alpha - c(X%*%beta))^2)/rchisq(1, J-1))
}
