% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Gaussian_Optim_Export.R
\name{DesignFit}
\alias{DesignFit}
\title{A Gaussian Process Functions}
\usage{
DesignFit(
  PatID,
  Data,
  Treat_M_List,
  Outcome_M_List = NULL,
  Outcome_Var,
  Treat_Vars,
  G_List,
  Psi,
  Normalized = TRUE,
  DR = FALSE,
  ...,
  Covtype,
  Numbr_Samp,
  IthetasU,
  IthetasL,
  IalphaU,
  IalphaL,
  Likelihood_Limits = NA,
  Prior_List = NULL,
  Prior_Der_List = NULL
)
}
\arguments{
\item{PatID}{patient identified}

\item{Data}{dataset used for analysis}

\item{Treat_M_List}{list of treatment models}

\item{Outcome_M_List}{list of outcome models}

\item{Outcome_Var}{outcome variable}

\item{G_List}{decision rule list}

\item{Psi}{points to explore in the parameter space}

\item{DR}{boolean indicating whether double robust estimator should be used}

\item{Covtype}{1 is for Matern 3/2; 2 is for Matern 5/2}

\item{Numbr_Samp}{number of random starts in the maximum likelihood}

\item{IthetasU}{upper limits in the bounded BFGS algorithm}

\item{IthetasL}{lower limits in the bounded BFGS algorithm}

\item{IalphaU}{lower limits in the bounded BFGS algorithm}

\item{IalphaL}{lower limits in the bounded BFGS algorithm}

\item{Likelihood_Limits}{limits for where to plot the parameters in the Gaussian Process Likelihood}

\item{Prior_List}{thsia}

\item{Prior_Der_List}{those}

\item{Treat_vars}{treatment variables}

\item{Bayes}{boolean indicating whether Bayesian methods will be used}

\item{Bayes_Seed}{seed to control random number generation}
}
\value{
This function returns the input parameters required for further fitting of the Gaussian Process when more samples are added

        as well as the current values of the parameters

This function allows you to fit a Bayesian MSM
}
\description{
This function fits Gaussian Process parameters to an initial set of design points.
}
