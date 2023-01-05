#-------------------------------------------------------------------------------------------------
#' Simulate data under a BMIM
#'
#' This function generates data under a BMIM under various assumptions. 
#' To be paired with main fitting function to investigate sensitivity to incorrect assumptions.
#'
#' @param N no. of observations
#' @param index_sizes M-vector of integers corresponding to index sizes. Corresponds to (L1,...L_M)
#' @param index_weights M-list of index weight vectors. (theta vectors of length Lm. will be transformed to have length 1) 
#' @param curve_funs M-list of index weight vectors. (theta vectors of length Lm)
#' @param interac_mat Upper triangular matrix with coefficients for interaction terms (i.e. 1,2 element of 0.5 means 0.5h1h2)
#' @param num_conf no. of (linear) confounders. Default 2
#' @param beta_conf vector of coefficients for confounders (length num_conf). Default 0s
#' @param corrs_within M-vector correlations (within index). Default 0
#' @param corrs_conf correlation among linear confounders. Default 0
#' @param sigma SD of error terms. Default 1
#' 
#' @importFrom mvtnorm
#' @return A list including samples of thetastar, theta, thetaPOS, rho, and posterior inclusion probabilities.
#' @author Glen McGee 
#' @export


simulate_bsmim2 <- function(N=50,
                            index_sizes=c(2,3),
                            index_weights=list(c(1,-1),
                                               c(3,2,1)),
                            curve_funs=list(function(x){x^2},
                                            function(x){sin(x)}),
                            interac_mat=NULL,
                            num_conf=2,
                            beta_conf=c(0,0),
                            corrs_within=0,
                            corrs_conf=0,
                            sigma=1){ 
  M <- length(index_sizes)
  if(length(index_weights)!=M){
    warning("incorrect weight length, using default")
    index_weights <- lapply(index_sizes,function(Lm){rep(1,Lm)})
  }
  if(length(curve_funs)!=M){
    warning("incorrect function length, using first")
    curve_funs <- lapply(curve_funs,function(fn){curve_funs[[1]]})
  }
  if(length(beta_conf)!=num_conf){
    warning("incorrect function length, using 0s")
    beta_conf <- rep(0,num_conf)
  }
  if(!is.null(interac_mat)){
    if(dim(interac_mat)[1]!=M |dim(interac_mat)[2]!=M){
      warning("incorrect interac_mat, using no interactions")
    }
  }
  
  ## L2 constraint
  theta <- lapply(index_weights,function(th){
    th/sqrt(sum(th^2))
  })
  
  X <- lapply(1:M, function(m){
    mvtnorm::rmvnorm(N,
                     rep(0,index_sizes[m]),
                     diag(1-corrs_within,index_sizes[m])+matrix(corrs_within,nrow=index_sizes[m],ncol=index_sizes[m]))
  })
  
  h_m <- lapply(1:M,function(m){
    curve_funs[[m]](as.numeric(X[[m]]%*%theta[[m]]))
  })
  
  ## main effects
  h <- rep(0,N)
  for(m in 1:M){
    h <- h+h_m[[m]]
  }
  ## add interactions
  if(!is.null(interac_mat)){ ## if exist
    for(m1 in 1:(M-1)){
      for(m2 in (m1+1):M){
        h <- h+interac_mat[m1,m2]*h_m[[m1]]*h_m[[m2]]
      }
    }
  }
  
  ## linear confounders
  z <- mvtnorm::rmvnorm(N,
                        rep(0,num_conf),
                        diag(1-corrs_conf,num_conf)+matrix(corrs_conf,nrow=num_conf,ncol=num_conf))

  ## generate outcomes
  y <- h+as.numeric(z%*%beta_conf)+rnorm(N,0,sigma)
  
  return(list(y=y,
              X=X,
              z=z))
}

