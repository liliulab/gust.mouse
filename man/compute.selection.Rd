% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r-functions.selection.R
\name{compute.selection}
\alias{compute.selection}
\title{Estimate selection coefficients of missense mutations and truncating mutations}
\usage{
compute.selection(obs_exp, joint = FALSE, indel = FALSE)
}
\arguments{
\item{obs_exp:}{A data frame with observed and expected counts}

\item{joint:}{A boolean flag to estimate joint probabilities}

\item{indel:}{A boolean flag to consider indel jointly with or independently from substitutions}
}
\value{
A vector with observed counts, end values and estimated probabilities
}
\description{
Estimate selection coefficients of missense mutations and truncating mutations
}
