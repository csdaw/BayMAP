bf <-
function(u_exp, u_mm, u_snp, q, dist, k, n, print.i = NULL){
  
  N <- length(n)
  l1 <- numeric(N)
  l2 <- numeric(N)
  
  if(dist == "binomial"){
    for(i in 1:N){
      l1[i] <- mean(unlist(lapply(u_exp, function(x){dbinom(k[i], n[i], x)})))
      l2[i] <- mean(unlist(apply(cbind(u_mm, u_snp, q), 1, function(x){x[3]*dbinom(k[i], n[i], x[1]) + (1-x[3])*dbinom(k[i], n[i], x[2])})))
      if(!is.null(print.i)){
        if((i - print.i) %% print.i == 0)
          print(i)
      }
      
    }
  }
  
  if(dist == "truncated"){
    for(i in 1:N){
      l1[i] <- mean(unlist(lapply(u_exp, function(x){dbinom(k[i], n[i], x)/(1-dbinom(0, n[i], x))})))
      l2[i] <- mean(unlist(apply(cbind(u_mm, u_snp, q), 1, function(x){x[3]*(dbinom(k[i], n[i], x[1])/(1-dbinom(0, n[i], x[1]))) + (1-x[3])*(dbinom(k[i], n[i], x[2])/(1-dbinom(0, n[i], x[2]))) })))
      if(!is.null(print.i)){
      if( (i - print.i) %%print.i == 0)
        print(i)
      }
    }
  }
  return(l1/l2)
}
