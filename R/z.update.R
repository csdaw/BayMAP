z.update <-
function(method, k, n, mu, p, q){
  
  w1 <- logpost_mu_single(mu=mu[1], method=method, k=k, n=n) + log(1-p) + log(q)
  w2 <- logpost_mu_single(mu=mu[2], method=method, k=k, n=n) + log(1-p) + log(1-q)
  w3 <- logpost_mu_single(mu=mu[3], method=method, k=k, n=n) + log(p)
  
  w.temp <- cbind(w1,w2,w3)
  w.temp <- exp(w.temp - apply(w.temp, 1, max))
  w <- w.temp/rowSums(w.temp)
  
  z.temp <- apply(w, 1, function(x)rmultinom(1,1,x))*c(1,2,3)
  colSums(z.temp)
  
}
