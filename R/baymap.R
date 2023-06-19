#' A Bayesian hierarchical model for the Analysis of PAR-CLIP data
#' 
#' @description The \code{baymap} function runs BayMAP on PAR-CLIP data to 
#'   detect PAR-CLIP induced T-to-C substitution positions on binding sites.
#' 
#' 
#' @param data a data frame with at least the count for mutations per genomic 
#'   position, the number of reads/coverage and the mutation type (e.g., T-to-C).
#' @param count the name of the variable that counts the number of mutations.
#' @param coverage the name of the variable that contains the number of reads.
#' @param mutation the name of the variable that contains the different types of mutations.
#' @param mutation.type the name of the mutation type that is induced by the PAR-CLIP method.
#' @param covariates a vector containing the names for the covariates for the 
#'   regression model, e.g., c("tpUTR", "cds", "fpUTR"). Intercept is automatically added as first variable.
#' @param dist the distribution for the number of mutations. Possible entries 
#'   are "truncated" (default) and "binomial.
#' @param dep a logical value for defining if dependencies between mismatches 
#'   and SNPs are considered (default) or not.
#' @param n.chains number of Markov chains (default: 1).
#' @param n.iter number of total iterations per chain (including burn in; 
#'   default 1500).
#' @param thin thinning rate. Must be a positive integer. Set \code{n.thin} > 1 
#'   to save memory and computation time if n.iter is large.
#' @param sd.mu a vector containing three values of standard deviations for the 
#'   sampling of \code{mu} with a normal jumping distribution.
#' @param inits.z a vector containing as inits an allocation for each position, 
#'   where 1 stands for an experimental induced substitution position, 2 for a 
#'   SNP and 3 for a mismatch.
#' @param inits.q a numerical value between 0 and 1 containing as init the 
#'   conditional probability for a mismatch position given the subsitions are 
#'   not experimentally induced. 
#' @param inits.mu a numerical vector containing as inits three values between 
#'   0 and 1 for the substitution probability due to the PAR-CLIP method, due to 
#'   SNPs and due to mismatches.
#' @param inits.beta a numerical vector containing as inits the parameter vector
#'   for the covariates. Only necessary if the vector \code{covariates} is not 
#'   \code{NULL}.
#' @param ran a logical value indicating if neighborhood dependencies should be 
#'   included via a random effect (default) or not.
#' @param cluster the name of the varialbe indicating to which cluster a 
#'   position belongs. Only necessary if \code{ran} is set to \code{TRUE}.
#' @param inits.tau a numerical value containing as inits the standard deviation
#'   of the random effect if \code{ran} is set to \code{TRUE}.
#' @param print.i a positive integer indicating if every ith iteration step 
#'   should be printed.
#' @param save_log a logical value indicating if temporary results should be 
#'   saved or not (default).
#' @param save_file file name where temporary results should be stored if 
#'   \code{save_log} is set to \code{TRUE}.
#'
#' @return The returned object is a list with sampled MCMC chains for each 
#'   parameter as entries and an entry with acceptance values for each sampled 
#'   value for the parameter \code{mu}. 
#' @references Huessler, E., Schaefer, M. Schwender, H., Landgraf, P. 
#'   (2019): BayMAP: A Bayesian hierarchical model for the analysis of PAR-CLIP 
#'   data. \emph{Bioinformatics}, 35(12), 1992-2000.
#' @author Eva-Maria Huessler, \email{eva-maria.huessler@uni-duesseldorf.de}
#' @seealso [predict.baymap()]
#' @export
#'
#' @examples
#' \dontrun{
#'   data(data_test)
#'   res <- baymap(data = data_test,
#'                 inits.mu = c(0.05, 0.85, 0.2), n.iter = 4500)
#' }
baymap <-
function(data, count = "count", coverage = "coverage", mutation = "mutation", 
                 mutation.type = "TC", covariates = NULL, dist = c("truncated", "binomial"), 
                 dep = TRUE, n.chains = 1, n.iter = 1500, thin = 1,
                 sd.mu = c(0.0001, 0.0001, 0.0001), inits.z = NULL, inits.q = NULL, 
                 inits.mu = NULL, inits.beta = NULL, ran = FALSE, cluster = "cluster", inits.tau = NULL,
                 print.i = NULL, save_log = FALSE, save_file = "./results_tmp.RData"){
    
            
            if(!is.data.frame(data)){
                stop("data must be a data frame.")
            }
            if(!is.character(count) | length(count) > 1){
                stop("The variable name of the count of mutations must be given as a character string of length 1.")
            }
            if(!is.character(coverage) | length(coverage) > 1){
                stop("The variable name of the read coverage must be given as a character string of length 1.")
            }
            if(!is.character(mutation) | length(mutation) > 1){
                stop("The variable name of the mutation type must be given as a character string of length 1.")
            }
            if(!is.character(mutation.type) | length(mutation.type) > 1){
                stop("The type of mutation that could be experimentally induced must be given as a character string of length 1.")
            }
            
            if(ran == FALSE){
                cluster = NULL
            }
    
            varnames <- c(count, coverage, mutation, covariates, cluster)
            varnames <- varnames[!is.na(varnames)]
            
            data <- data[ , varnames]
            
            k <- data[,count]
            n <- data[,coverage]
            TC <- data[,mutation] == mutation.type
            if(ran == FALSE){
                if(is.null(covariates)){
                    X <- as.matrix(rep(1, sum(TC)))
                }else{
                    X <- as.matrix(data.frame(1, data[TC, covariates]))
                }
            }else{
                zone <- data[TC, cluster]
                if(is.null(covariates)){
                    X <- matrix(numeric(length(unique(data[TC, cluster])) * 1), ncol = 1)
                }else{
                    X <- matrix(numeric(length(unique(data[TC, cluster])) * (length(covariates)+1)), ncol = (length(covariates)+1))
                }
                l <- 1
                for(i in unique(data[TC, cluster])){
                    X[l,] <- cbind(1, data[(which(data[, cluster] == i) & TC), covariates][1])
                    l <- l+1
                }
            }
            dist <- match.arg(dist)
            if(dep){
                bm <- "_bm"
            }else{
                bm <- ""
            }
            if(dist == "truncated"){
                dist <- "zero"
            }
            method <- paste0(dist, bm)
            
       
        sims <- array(NA, c(n.iter, n.chains))
        z.post <- array(NA, c(n.iter, length(k), n.chains))
        mu.post <- array(NA, c(n.iter, 3, n.chains))
        q.post <- array(NA, c(n.iter, n.chains))
        p <- numeric(length(k))
        if(ran){
            nzone <- table(zone)[rank(unique(zone))]
            J <-length(nzone)
            alpha.post <- array(NA, c(n.iter, J, n.chains))
            beta.post <- array(NA, c(n.iter, ncol(X),n.chains))
            tau.post <- array(NA, c(n.iter, n.chains))
        }else{
            beta.post <- array(NA, c(n.iter, dim(X)[2], n.chains))
            beta <- numeric(dim(X)[2])
        }
        acc.post <- list(mu = array(NA, c(n.iter*thin, 3, n.chains)))
        
        for(m in 1:n.chains){
            if(ran){
                beta <- rnorm(ncol(X), 0, 1)
                tau <- 1
                alpha <- numeric(J)
                p[TC == 1] <- pnorm(rep(alpha, nzone), sd = 1)
            }else{
                if(!is.null(inits.beta)){
                    beta <- inits.beta
                }else{beta <- rnorm(dim(X)[2],0,1)}
                p[TC ==1] <- pnorm(X%*%beta, sd=1)
            }
            
            if(!is.null(inits.q)){
                q <- inits.q
            }else{q <- rbeta(1,1,1)}
            
            rq <- cbind((1-p)*q, (1-p)*(1-q),p)
            
            if(!is.null(inits.z)){
                z <- inits.z
            }else{
                z.temp <- apply(rq, 1, function(x)rmultinom(1,1,x))*c(1,2,3)
                z <- colSums(z.temp)
            }
    
            if(!is.null(inits.mu)){
                mu <- inits.mu
            }else{
                mu <- numeric(3)
                mu[1] <- qbeta(runif(1, min=pbeta(0, 1, 1), max=pbeta(0.25, 1 , 1)), 1, 1)
                mu[2] <- 1-3*mu[1]
                mu[3] <- rbeta(1,1,1)
            }
            
    
            if(method != "binomial"){cur_logpost_mu <- logpost_mu(mu=mu, method=method, k=k, n=n, z=z)}
            
            counter <- 1
            
            for(t in 1:(n.iter*thin)){
                mu.res <- mu.update(method=method, k=k, n=n, z=z, mu=mu, cur_logpost_mu=cur_logpost_mu, sd.mu = sd.mu)
                mu <- mu.res$est
                z <- z.update(method=method, k=k, n=n, mu=mu, p=p, q=q)
                q <- q.update(z=z)
                if(ran){
                    y <- y.update.ran(alpha = alpha, nzone = nzone, z=z, TC=TC)
                    alpha <- alpha.update(y = y, zone = zone, nzone = nzone, J = J, beta = beta, X=X, tau=tau)
                    beta <- beta.update.ran(alpha = alpha, tau = tau, X=X)
                    tau <- tau.update(alpha, beta, X=X, J)
                }else{
                    y <- y.update(X=X, beta=beta, z=z, TC=TC)
                    beta.res <- beta.update(beta=beta, X=X, z=z, y=y)
                    beta <- beta.res$est
                }
    
                if(ran){
                    p[TC == 1] <- pnorm(rep(alpha, nzone), sd=1)
                }else{
                    p[TC ==1] <- pnorm(X%*%beta, sd=1)
                }
                
                if(!is.null(print.i)){   
                    if((t - print.i) %% (print.i*thin) == 0){
                        print(paste0("chain number ",m, ", sim number ", counter))  }
                }
                if((t-1) %% thin == 0){
                    mu.post[counter, ,m] <- mu
                    z.post[counter,,m] <- z
                    q.post[counter,m] <- q
                    if(ran){
                        alpha.post[counter, ,m] <- alpha-X%*%beta
                        beta.post[counter, ,m] <- beta
                        tau.post[counter,m] <- tau
                    }else{
                        beta.post[counter, ,m] <- beta
                    }
    
                    counter <- counter + 1
                    
                }
                if((t-1) %% (thin*500) == 0){
                    if(save_log == TRUE){
                        if(ran){
                             results_tmp <- list("mu" = mu.post, "z" = z.post, "q" = q.post, "alpha" = alpha.post, "beta" = beta.post, "tau" = tau.post)
                        }else{
                            results_tmp <- list("mu" = mu.post, "z" = z.post, "q" = q.post, "beta" = beta.post)
                        }
                        save(results_tmp, file = save_file)
                    }
                }
                
                if(ran == FALSE){
                    acc.post$beta[t, ,m] <- beta.res$acc
                }
                acc.post$mu[t, , m] <- mu.res$acc
            }
        }
        if(ran){
            results <- list("mu" = mu.post, "z" = z.post, "q" = q.post, "alpha" = alpha.post, "beta" = beta.post, "tau" = tau.post, "acceptance" = acc.post)
        }else{
            results <- list("mu" = mu.post, "z" = z.post, "q" = q.post, "beta" = beta.post, "acceptance" = acc.post)
        }
        class(results) <- "baymap"
        return(results)
    }
