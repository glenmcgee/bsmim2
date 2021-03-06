% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bsmim_predict2.R
\name{predict_hnew_assoc2}
\alias{predict_hnew_assoc2}
\title{Compute exposure-associations: hnew(x*)-hnew(median(x)), setting all other exposures to median}
\usage{
predict_hnew_assoc2(
  object,
  qtl_lims = c(0.05, 0.95),
  qtls = list(NA),
  compare_qtl = 0.5,
  overall = FALSE,
  approx = FALSE
)
}
\arguments{
\item{object}{An object of class bsmim}

\item{qtl_lims}{Quantile limits for grid of exposure values}

\item{qtls}{M-list of Lm-vectors containing quantiles to set exposure components to}

\item{compare_qtl}{quantile to compare all others to}

\item{overall}{should the 'overall' effect be reported}

\item{approx}{should the approximate method be used (not recommended)}
}
\value{
a dataframe containing predicted values.
}
\description{
Compute exposure-associations: hnew(x*)-hnew(median(x)), setting all other exposures to median
}
