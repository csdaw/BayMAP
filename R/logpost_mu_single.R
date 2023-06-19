logpost_mu_single <-
function(mu, method, k, n){
  like <- k*log(mu) + (n-k)*log(1-mu) 
  if(method == "zero_bm" || method == "zero"){
    like <- like - log(1 - (1-mu)^n)
  }
  return(like)
}
