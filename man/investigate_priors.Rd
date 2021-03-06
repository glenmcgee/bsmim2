% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/investigate_priors.R
\name{investigate_priors}
\alias{investigate_priors}
\title{Investigate priors on thetastar}
\usage{
investigate_priors(
  num_components = 8,
  R = 2000,
  constraint = 0,
  spike_slab = TRUE,
  gauss_prior = TRUE,
  prior_theta_slab_sd = 0.25,
  prior_theta_slab_bounds = c(0, 100),
  prior_pi = c(1, 1),
  prior_slabpos = c(0.4, 1.6),
  prior_slabpos_shape_inf = NULL,
  prior_alphas = NULL,
  prior_slabrho = c(4, 2),
  horseshoe = FALSE,
  kappa = 1,
  prior_tau = 0.01
)
}
\arguments{
\item{num_components}{no. components in index}

\item{R}{no. of samples}

\item{constraint}{Which constraints should it use? 0=unconstrained, 1= non-negative, 2=positive+dirichlet}

\item{spike_slab}{Use spike and slab for variable selection (new version assumes this is true)}

\item{gauss_prior}{gaussian priors on thetastar ## only if spike_slab=TRUE}

\item{prior_theta_slab_sd}{sd of slab}

\item{prior_theta_slab_bounds}{bounds of unif for invunif prior on thetastar^2 (for gauss_prior=FALSE)}

\item{prior_pi}{prior for pi in spike & slab}

\item{prior_slabpos}{shape and rate for gamma prior on thetastar (under constraints=1)}

\item{prior_slabpos_shape_inf}{vector of multiple shape params for gamma prior on thetastar (under constraints=1) when using informative priors with variable selection. Overrides shape param in prior_slabpos[1]}

\item{prior_alphas}{vector of alpha hyperparameters for dirichlet prior on weights (under constraints=2)}

\item{prior_slabrho}{shape and rate for gamma prior on rho^{1/2} (under constraints=2)}

\item{horseshoe}{Old version with spike_slab=FALSE: Should we use the horseshoe prior for indexwise selection (default FALSE)}

\item{kappa}{Old version with spike_slab=FALSE: scale parameter kappa used in non-horseshoe (sd of thetastar is kappa * nu, nu is chi^2_1)}

\item{prior_tau}{Old version with spike_slab=FALSE: Hyperparameter \eqn{\tau_{0}^2}, scale parameter for half Cauchy prior for \eqn{\tau}. (sd of thetastar is tau* nu, where tau is half-cauchy(0,tau_0^2) and nu is half-cauchy(0,1))}
}
\value{
A list including samples of thetastar, theta, thetaPOS, rho, and posterior inclusion probabilities.
}
\description{
This function samples from the priors for thetastar/theta/thetaPOS (different parameterizations) to facilitate choosing priors.
Note: only considers a weights in a single index. As such, arguments that were M-lists of Lm vectors in bsmim2() are now just Lm vectors etc
}
\author{
Glen McGee
}
