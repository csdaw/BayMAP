#' Predict method for BayMAP results
#' 
#' @description Predictions if T-to-C substitution positions are PAR-CLIP 
#'   induced substitutions or not. Results of several PAR-CLIP experiments can be 
#'   combined for a combined prediction.
#'
#' @param object either a \code{baymap} object or a list of \code{baymap} 
#'   objects if several PAR-CLIP experiments should be analyzed jointly. If a list
#'   with \code{baymap} objects is read in, the class "baymap" should be assigned
#'   to the list prior to the analysis by \code{\link[base]{class}}.
#' @param data either a data frame with at least the count for mutations per 
#'   genomic position, the number of reads/coverage and the mutation type 
#'   (e.g., T-to-C) or a list of data frames.
#' @param chr the name of the variable that contains the chromosome information.
#' @param pos the name of the variable that contains the position information.
#' @param count the name of the variable that counts the number of mutations.
#' @param coverage the name of the variable that contains the number of reads.
#' @param mutation the name of the variable that contains the different types of mutations.
#' @param mutation.type the name of the mutation type that is induced by the PAR-CLIP method.
#' @param covariates a vector containing the names for the covariates for the 
#'   regression model, e.g., c("tpUTR", "cds", "fpUTR"). Intercept is 
#'   automatically added as first variable.
#' @param dist the distribution for the number of mutations. Possible 
#'   entries are "truncated" (default) and "binomial.
#' @param ran a logical value indicating if neighborhood dependencies should be 
#'   included via a random effect (default) or not.
#' @param cluster the name of the varialbe indicating to which cluster a 
#'   position belongs. Only necessary if \code{ran} is set to \code{TRUE}.
#' @param print.i a positive integer indicating if every ith iteration step 
#'   should be printed.
#' @param thin an additional thinning rate that should be applied on the 
#'   \code{baymap} outcome. Must be a positive integer.
#' @param burn.in length of burn in, i.e. number of iterations to discard at the
#'   beginning.
#' @param ... further arguments for \code{predict}.
#'
#' @return If a single PAR-CLIP data set is analyzed, the returned object is a 
#'   data frame including the prior odds, the bayes factor and the posterior odds.
#'   If several PAR-CLIP data sets are analyzed, the returned object is a list 
#'   combined predictions as well as separate predictions for each data set. 
#'   Posterior odds greater than one means that the probability of having a 
#'   method-induced substitution position given the data is larger than 0.5.
#' 
#' @references Huessler, E., Schaefer, M. Schwender, H., Landgraf, P. (2019): 
#'   BayMAP: A Bayesian hierarchical model for the analysis of PAR-CLIP data. 
#'   \emph{Bioinformatics}, 35(12), 1992-2000.
#' @author Eva-Maria Huessler, \email{eva-maria.huessler@uni-duesseldorf.de}
#' @note Predictions of the separate analyses are made for all entries of the 
#'   included data set even for substitution types other than T-to-C.
#' @seealso [baymap()]
#'   
#' @export
#'
#' @examples
#' \dontrun{
#'   data(data_test)
#'   res <- baymap(data = data_test,
#'                 inits.mu = c(0.05, 0.85, 0.2), n.iter = 4500)
#'   data_new <- predict(res, data_test, burn.in = 3000)
#' }
predict.baymap <-
  function(object, data, chr = "chr", pos = "pos", count = "count", coverage = "coverage", mutation = "mutation", 
           mutation.type = "TC", covariates = NULL, dist = c("truncated", "binomial"),
           ran = FALSE, cluster = NULL, print.i = 100, thin = NULL, burn.in = 0, ...){
    
    if(is.data.frame(data)){
      data <- list(data)
      object <- list(object)
    }
    
    if(!is.null(cluster)){
      cluster.id <- cluster
    }
    
    if(!is.character(count) | length(count) > 1){
      stop("The variable name of the count of mutations must be given as a character string of length 1.")
    }
    if(!is.character(coverage) | length(coverage) > 1){
      stop("The variable name of the read coverage must be given as a character string of length 1.")
    }
    # Number of data sets
    n <- length(data)
    
    if(length(object) != n){
      stop("Object and data must have the same length.")
    }
    
    dist <- match.arg(dist)  
    
    if(ran == FALSE){
      alpha = NULL
    }else{
      alpha <- vector("list", n)
    }
    
    if(is.null(thin)){
      thin <- 1
    }
    
    for(i in 1:n){
      t <- seq(burn.in + 1, nrow(object[[i]]$mu), by = thin)
      if(ran){
        alpha[[i]] <- object[[i]]$alpha[t,,1]
        
        alpha_rep <- matrix(numeric(length(t)*nrow(data[[i]])), nrow = length(t))
        cluster <- unique(data[[i]][data[[i]][,mutation] == mutation.type,cluster.id])
        for(j in 1:nrow(data[[i]][data[[i]][,mutation] == mutation.type,])){
          cl <- data[[i]][data[[i]][,mutation] == mutation.type,]$cluster[j][]
          if(cl != 0){
            alpha_rep[,j] <- alpha[[i]][, which(cluster == cl)]
          }
        }
      }else{
        alpha_rep <- NULL
      }
      
      data[[i]]$BayesFactor <- bf(u_exp = object[[i]]$mu[t,3,1], u_mm  = object[[i]]$mu[t,1,1], u_snp = object[[i]]$mu[t,2,1], q  = object[[i]]$q[t,1], k = data[[i]][,count], n = data[[i]][,coverage], dist = dist, print.i = print.i)
      if(!is.null(covariates)){
        X <- cbind(1, data[[i]][, covariates])
      }else{
        X <- matrix(1, nrow = nrow(data[[i]]))
      }
      data[[i]]$PriorOdds <- prO(X = X, beta = object[[i]]$beta[t,,1], alpha = alpha_rep)
      data[[i]]$PostOdds <- data[[i]]$BayesFactor*data[[i]]$PriorOdds
    }
    
    if(n == 1){
      return(data[[1]])
    }else{
      data_TC_tmp <- vector("list", n)
      for(i in 1:n){
        if(!is.null(covariates)){
          colnames(data[[i]])[!(colnames(data[[i]]) %in% c(chr, pos, mutation, covariates))] <- paste(colnames(data[[i]])[!(colnames(data[[i]]) %in% c(chr, pos, mutation, covariates)) ], i, sep = "_")
        }else{
          colnames(data[[i]])[!(colnames(data[[i]]) %in% c(chr, pos, mutation)) ] <- paste(colnames(data[[i]])[!(colnames(data[[i]]) %in% c(chr, mutation, pos)) ], i, sep = "_")
        }
        data_TC_tmp[[i]] <- data[[i]][data[[i]][, mutation] == mutation.type, ]
      }
      if(!is.null(covariates)){
        data_TC_merge <- merge(all = TRUE, suffixes = c("", ""), by = c(chr, pos, mutation, covariates), data_TC_tmp[[1]] , data_TC_tmp[[2]])
      }else{
        data_TC_merge <- merge(all = TRUE, suffixes = c("", ""), by = c(chr, pos), data_TC_tmp[[1]] , data_TC_tmp[[2]])
      }
      if(n >= 3){
        if(!is.null(covariates)){
          for(i in 1:(n-2)){ 
            data_TC_merge <- merge(all = TRUE, suffixes = c("", ""), by = c(chr, pos, covariates), data_TC_merge , data_TC_tmp[[i+2]])
          }
        }else{
          for(i in 1:(n-2)){
            data_TC_merge <- merge(all = TRUE, suffixes = c("", ""), by = c(chr, pos, covariates), data_TC_merge , data_TC_tmp[[i+2]])
          }
        }
      }
      cluster_new <- vector("list", n)
      cl_ges <- vector("list", n)
      alpha_rep <- vector("list", n)
      beta <- vector("list", n)
      X <- cbind(1, data_TC_merge[, covariates])
      for(j in 1:n){
        t <- seq(burn.in + 1, nrow(object[[i]]$mu), by = thin)
        m <- matrix(rep(NA, length(t)*nrow(data_TC_merge)), nrow = length(t))
        alpha_rep[[j]] <- m
        beta[[j]] <- object[[j]]$beta[t,,1]
        alpha[[j]] <- object[[j]]$alpha[t,,1]
        cluster_new[[j]] <- unique(data_TC_tmp[[j]][,paste(cluster.id, j, sep = "_")])
        cl_ges[[j]] <- data_TC_merge[,paste0(cluster.id, "_", j)]
      }
      for(l in 1:nrow(data_TC_merge)){
        for(k in 1:length(alpha_rep)){
          cl <- cl_ges[[k]][l]
          if(!is.na(cl)){
            alpha_rep[[k]][,l] <- alpha[[k]][, which(cluster_new[[k]] == cl)]
          }
        }
      }
      
      data_TC_merge$PriorOdds_combined <- prO(X = X, beta = beta, alpha = alpha_rep)
      data_TC_merge$PostOdds_combined <- data_TC_merge$BayesFactor_1^(!is.na(data_TC_merge$BayesFactor_1)) * data_TC_merge$BayesFactor_2^(!is.na(data_TC_merge$BayesFactor_2))
      if(n >= 3){
        for(i in (n-2)){
          bf_tmp <- eval(parse(text = paste0("data_TC_merge$BayesFactor_", (i+2))))
          data_TC_merge$PostOdds_combined <- data_TC_merge$PostOdds_combined * bf_tmp^(!is.na(bf_tmp)) 
        }
      }
      data_TC_merge$PostOdds_combined <- data_TC_merge$PostOdds_combined * data_TC_merge$PriorOdds_combined
      
      return(list(single_data = data, combined_data = data_TC_merge))
    }
  }
