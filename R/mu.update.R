mu.update <-
function(method, k, n, z, mu, cur_logpost_mu, sd.mu){
  acc <- numeric(3)
  
  if(method == "binomial"){
    mu[1]  <- rbeta(1, sum(k[z==1]) + 1, sum(n[z==1]) - sum(k[z==1]) +1)
    mu[2] <- rbeta(1, sum(k[z==2]) + 1, sum(n[z==2]) - sum(k[z==2]) +1)
    mu[3] <- rbeta(1, sum(k[z==3]) + 1, sum(n[z==3]) - sum(k[z==3]) +1)
  }
  
  if(method == "binomial_bm"){
    cur_logpost_mu <- logpost_mu(mu=mu, method=method, k=k, n=n, z=z)
    mu.star <- mu
    mu.star[1] <- rnorm(n=1, mean=mu[1], sd=sd.mu[1])
    mu.star[2] <- 1-3*mu.star[1]
    can_logpost_mu <- logpost_mu(mu=mu.star, method=method, k=k, n=n, z=z)
    R <- min(1, exp(can_logpost_mu - cur_logpost_mu))
    U <- runif(1)
    if(U<R){
      mu <- mu.star
      cur_logpost_mu <- can_logpost_mu
      acc[1] <- 1
    }
    
    mu[3] <- rbeta(1, sum(k[z==3]) + 1, sum(n[z==3]) - sum(k[z==3]) +1)
  }
  
  if(method == "zero"){
    cur_logpost_mu <- logpost_mu(mu=mu, method=method, k=k, n=n, z=z)
    for(i in 1:3){
      mu.star <- mu
      mu.star[i] <-rnorm(n=1, mean = mu[i], sd=sd.mu[i])
      can_logpost_mu <- logpost_mu(mu=mu.star, method=method, k=k, n=n, z=z)
      R <- min(1, exp(can_logpost_mu - cur_logpost_mu))
      U <- runif(1)
      if(U<R){
        mu <- mu.star
        cur_logpost_mu <- can_logpost_mu
        acc[i] <- 1
      }
    }
  }
  
  if(method == "zero_bm"){
    cur_logpost_mu <- logpost_mu(mu=mu, method=method, k=k, n=n, z=z)
    mu.star <- mu
    mu.star[1] <- rnorm(n=1, mean=mu[1], sd=sd.mu[1])
    mu.star[2] <- 1-3*mu.star[1]
    can_logpost_mu <- logpost_mu(mu=mu.star, method=method, k=k, n=n, z=z)
    R <- min(1, exp(can_logpost_mu - cur_logpost_mu))
    U <- runif(1)
    if(U<R){
      mu <- mu.star
      cur_logpost_mu <- can_logpost_mu
      acc[1] <- 1
      acc[2] <- 1
    }
    mu.star <- mu
    mu.star[3] <- rnorm(n=1, mean=mu[3], sd=sd.mu[3])
    can_logpost_mu <- logpost_mu(mu=mu.star, method=method, k=k, n=n, z=z)
    R <- min(1, exp(can_logpost_mu - cur_logpost_mu))
    U <- runif(1)
    if(U<R){
      mu <- mu.star
      cur_logpost_mu <- can_logpost_mu
      acc[3] <- 1
    }
  }
  l <- list(est = mu, acc = acc)
  return(l)
}
