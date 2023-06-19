logpost_mu <-
function(mu, method, k, n, z){
  if(method == "binomial_bm"){
    if(!any(mu < 0)){
      like <- sum(lchoose(n, k) 
                  + sum(k[z==1])*log(mu[1]) + sum(k[z==2])*log(mu[2]) + sum(k[z==3])*log(mu[3])
                  + (sum(n[z==1]) - sum(k[z==1])) * log(1 - mu[1]) + (sum(n[z==2]) - sum(k[z==2])) * log(1 - mu[2]) + (sum(n[z==3]) - sum(k[z==3])) * log(mu[3])
      )
    }else{like <- -Inf}
    prior <- dunif(mu[1], min= 0, max = 0.25, log=TRUE) +  dunif(mu[2], min= 0.25, max = 1, log=TRUE) + dunif(mu[3], log=TRUE)
  }
  
  if(method == "zero"){
    if(!any(mu < 0)){
      like <- sum(lchoose(n, k) 
                  + sum(k[z==1])*log(mu[1]) + sum(k[z==2])*log(mu[2]) + sum(k[z==3])*log(mu[3])
                  + (sum(n[z==1]) - sum(k[z==1])) * log(1 - mu[1]) + (sum(n[z==2]) - sum(k[z==2])) * log(1 - mu[2]) + (sum(n[z==3]) - sum(k[z==3])) * log(1 - mu[3])
                  - sum(log(1-(1-mu[z])^n)))
    }else{like <- -Inf}
    prior <- dunif(mu[1], log=TRUE) + dunif(mu[2], log=TRUE) + dunif(mu[3], log=TRUE)
  }
  
  if(method == "zero_bm"){
    if(!any(mu < 0)){
      like <- sum(lchoose(n, k) 
                  + sum(k[z==1])*log(mu[1]) + sum(k[z==2])*log(mu[2]) + sum(k[z==3])*log(mu[3])
                  + (sum(n[z==1]) - sum(k[z==1])) * log(1 - mu[1]) + (sum(n[z==2]) - sum(k[z==2])) * log(1 - mu[2]) + (sum(n[z==3]) - sum(k[z==3])) * log(1 - mu[3])
                  - sum(log(1-(1-mu[z])^n)))
    }else{like <- -Inf}
    prior <- dunif(mu[1], min= 0, max = 0.25, log=TRUE) +  dunif(mu[2], min= 0.25, max = 1, log=TRUE) + dunif(mu[3], log=TRUE)
  }
  
  return(like +prior)
}
