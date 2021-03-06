% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BivariateAutoregressiveHMM.R
\name{CalcTran}
\alias{CalcTran}
\title{Calculates latent state transition probabilities}
\usage{
CalcTran(
  forw,
  backw,
  tran,
  algae_data,
  toxin_data,
  mu_a,
  k_a,
  tox_em,
  states,
  toxin_states
)
}
\arguments{
\item{forw}{Forward probabilities}

\item{backw}{Backward probabilities}

\item{tran}{Old Markov state transition probabilities}

\item{algae_data}{Algae count data}

\item{toxin_data}{Toxin data measured in micro grams / OA equivalent}

\item{mu_a}{mean negative binomial parameter}

\item{k_a}{size negative binomial parameter}

\item{tox_em}{Classification matrix from ordinal regression parameters}
}
\value{
Returns initial probabilities
}
\description{
Calculates latent state transition probabilities
}
