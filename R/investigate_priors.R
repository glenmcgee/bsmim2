
#-------------------------------------------------------------------------------------------------
#' Investigate priors on thetastar
#'
#' This function samples from the priors for thetastar/theta/thetaPOS (different parameterizations) to facilitate choosing priors.
#' Note: only considers a weights in a single index. As such, arguments that were M-lists of Lm vectors in bsmim2() are now just Lm vectors etc
#'
#' @param num_components no. components in index
#' @param R no. of samples
#' @param constraint Which constraints should it use? 0=unconstrained, 1= non-negative, 2=positive+dirichlet
#' @param spike_slab Use spike and slab for variable selection (new version assumes this is true)
#' @param gauss_prior gaussian priors on thetastar ## only if spike_slab=TRUE
#' @param prior_theta_slab_sd sd of slab
#' @param prior_theta_slab_bounds bounds of unif for invunif prior on thetastar^2 (for gauss_prior=FALSE)
#' @param prior_pi prior for pi in spike & slab
#' @param prior_slabpos shape and rate for gamma prior on thetastar (under constraints=1)
#' @param prior_slabpos_shape_inf vector of multiple shape params for gamma prior on thetastar (under constraints=1) when using informative priors with variable selection. Overrides shape param in prior_slabpos[1]
#' @param prior_alphas vector of alpha hyperparameters for dirichlet prior on weights (under constraints=2)
#' @param prior_slabrho shape and rate for gamma prior on rho^{1/2} (under constraints=2)
#' @param horseshoe Old version with spike_slab=FALSE: Should we use the horseshoe prior for indexwise selection (default FALSE)
#' @param kappa Old version with spike_slab=FALSE: scale parameter kappa used in non-horseshoe (sd of thetastar is kappa * nu, nu is chi^2_1)
#' @param prior_tau Old version with spike_slab=FALSE: Hyperparameter \eqn{\tau_{0}^2}, scale parameter for half Cauchy prior for \eqn{\tau}. (sd of thetastar is tau* nu, where tau is half-cauchy(0,tau_0^2) and nu is half-cauchy(0,1))
#' @importFrom stats rbeta rbinom rnorm runif rcauchy rchisq rgamma
#' @return A list including samples of thetastar, theta, thetaPOS, rho, and posterior inclusion probabilities.
#' @author Glen McGee 
#' @export


investigate_priors <- function(num_components=8, # no. components in index
                               R=2000, ## no. of samples
                               constraint=0, ## 0=unconstrained, 1= non-negative, 2=positive+dirichlet
                               spike_slab=TRUE, ## spike and slab for variable selection (new version assumes this is true)
                               gauss_prior=TRUE, ## gaussian priors on thetastar ## only if spike_slab=TRUE
                               prior_theta_slab_sd=0.25, ## sd of slab
                               prior_theta_slab_bounds=c(0,100), ## bounds of unif for invunif prior on thetastar^2 (for gauss_prior=FALSE)
                               prior_pi=c(1,1), ## prior for pi in spike & slab
                               prior_slabpos=c(0.4,1.6), ## shape and rate for gamma prior on thetastar (under constraints=1)
                               prior_slabpos_shape_inf=NULL, ## vector of multiple shape params for gamma prior on thetastar (under constraints=1) when using informative priors with variable selection. Overrides shape param in prior_slabpos[1]
                               prior_alphas=NULL,    ## vector of alpha hyperparameters for dirichlet prior on weights (under constraints=2)
                               prior_slabrho=c(4,2), ## shape and rate for gamma prior on rho^{1/2} (under constraints=2)
                               ## OLD VERSION with no spike and slab prior:
                               horseshoe=FALSE, ## Use the horseshoe prior for indexwise selection (default FALSE)
                               kappa=1, ## scale parameter kappa used in non-horseshoe (sd of thetastar is kappa * nu, nu is chi^2_1)
                               prior_tau=0.01 ){ ## Hyperparameter \eqn{\tau_{0}^2}, scale parameter for half Cauchy prior for \eqn{\tau}. (sd of thetastar is tau* nu, where tau is half-cauchy(0,tau_0^2) and nu is half-cauchy(0,1))
  
  
  
  thetastar <- matrix(NA,nrow=R,ncol=num_components)  
  if(constraint==0){ ## unconstrained index
    
    if(spike_slab==TRUE){
      
      
      if(gauss_prior==TRUE){ ## gaussian prior on thetastar ( this is standard in the partially constrained version)
        for(l in 1:num_components){
          ## for variable selection
          pi <- rbeta(R,prior_pi[1],prior_pi[2])
          delta <- rbinom(R,1,pi)
          
          slab <- rnorm(R,0,prior_theta_slab_sd)
          thetastar[,l] <- slab*delta+0
          
        }
      }else{ ## else inverse uniform prior on thetastar^2 (only implemented when not using any constraints)
        for(l in 1:num_components){
          ## for variable selection
          pi <- rbeta(R,prior_pi[1],prior_pi[2])
          delta <- rbinom(R,1,pi)
          
          slab <- sqrt(1/runif(R,prior_theta_slab_bounds[1],prior_theta_slab_bounds[2]))
          flip <- rbinom(R,1,0.5) ## half prob of being negative (since it is symmetric)
          slab[flip==1] <- -slab[flip==1]
          thetastar[,l] <- slab*delta+0 ## var selection
          
        }
      }  
    }else{ ## no spike and slab (old version)
      
      if(horseshoe==TRUE){
        tau <- abs(rcauchy(R,0,prior_tau)) # half cauchy
        for(l in 1:num_components){
          
          nu <- abs(rcauchy(R,0,1)) # half cauchy
          thetastar[,l] <- rnorm(R,0,tau*nu) ## variance of nu_{ml}^2 * tau^2
          
        }
      }else{ ## no horseshoe 
        nu <- rchisq(R,1) # chi-square 1 for nu. same nu for full index
        for(l in 1:num_components){
          
          thetastar[,l] <- rnorm(R,0,kappa*nu) ## variance of kappa^2 * nu_m^2
          
        }
      }
    }
    
    
    ## original scale thetas
    rho <- apply(thetastar,1,function(x) sum(x^2))
    theta <- thetastar/sqrt(rho)
    
    ## handle all 0s
    theta[is.na(theta)] <- 0 
    
    ## identifiability constraint (sum is non-negative)
    cmw <- rowMeans(theta) 
    if(any(cmw<0)){
      theta[which(cmw<0),] <- -theta[which(cmw<0),]
    }
    
    ## thetaPOS not defined
    thetaPOS <- NULL
    
    
    
  }else if(constraint==1){ # always use spike and slab
    
    
    if(is.null(prior_slabpos_shape_inf)){
      prior_slabpos_shape_inf <- rep(prior_slabpos[1],num_components)
    }else if(length(prior_slabpos_shape_inf)!=num_components){
      print(paste("prior_slabpos_shape_inf should be an Lm-vector for this function"))
      prior_slabpos_shape_inf <- rep(prior_slabpos[1],num_components)
    }
    
    for(l in 1:num_components){
      pi <- rbeta(R,prior_pi[1],prior_pi[2])
      delta <- rbinom(R,1,pi)
      slab <- rgamma(R,shape=prior_slabpos_shape_inf[l],rate=prior_slabpos[2])
      thetastar[,l] <- slab*delta+0
    }
    
    ## original scale thetas
    rho <- apply(thetastar,1,function(x) sum(x^2))
    theta <- thetastar/sqrt(rho)
    
    ## handle all 0s
    theta[is.na(theta)] <- 0 
    
    ## thetaPOS--proportion scale
    thetaPOS <- thetastar/apply(thetastar,1,sum)
    thetaPOS[is.na(thetaPOS)] <- 0 
    thetaPOS[is.infinite(thetaPOS)] <- 0 
    
  }else if(constraint==2){
    if(is.null(prior_alphas)){
      prior_alphas <- rep(1,num_components)
    }else if(length(prior_alphas)!=num_components){
      print(paste("prior_alphas should be an Lm-vector for this function"))
      prior_alphas <- rep(1,num_components)
    }
    
    x <- matrix(NA,nrow=R,ncol=num_components)
    ## dirichlet from exponential
    for(l in 1:num_components){
      x[,l] <- rgamma(R,prior_alphas[l],1)
    }
    thetaPOS <- x/apply(x,1,sum) ## standardize for dirichlet
    
    ## separately draw \rho'^{1/2}
    sqrt_rhoprime <- rgamma(R,shape=prior_slabrho[1],rate=prior_slabrho[2])
    
    ## reconstruct thetastar
    thetastar <- sqrt_rhoprime*thetaPOS
    
    ## original scale thetas
    rho <- apply(thetastar,1,function(x) sum(x^2))
    theta <- thetastar/sqrt(rho)
    
  }
  
  ## compute posterior inclusion probabilities
  PIP <- apply(thetastar,2,function(x)100*mean(x!=0))
  
  return(list(thetastar=thetastar,
              theta=theta,
              rho=rho,
              thetaPOS=thetaPOS,
              PIP=PIP))
  
}


