prO <-
function(X, beta, alpha = NULL){
  if(is.matrix(beta)){
    beta <- list(beta)
  }
  if(is.matrix(alpha)){
    alpha <- list(alpha)
  }
  xb <- lapply(beta, function(x){x %*% t(X)})
  if(!is.null(alpha)){
    for(j in 1:length(alpha)){
      xb[[j]] <- xb[[j]] + alpha[[j]]
    }
  }
  pi <- sapply(xb, function(x){colMeans(pnorm(x, 0, 1))})
  pi <- rowMeans(pi)
  odds <- pi/(1-pi)
  return(odds)
}
