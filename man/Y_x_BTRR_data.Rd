% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/001_Y_x_BTRR_data.R
\name{Y_x_BTRR_data}
\alias{Y_x_BTRR_data}
\title{Bayesian Tensor Response Regression Data Generation}
\usage{
Y_x_BTRR_data(tens_dim, subjects, n.time, margin_size, CNR, obs.var)
}
\arguments{
\item{tens_dim}{A scalar (e.g. 2 or 3) for the dimension of the generated response data for a given time and subject}

\item{subjects}{A scalar with the sample size of the generated data}

\item{n.time}{A scalar with the number of time steps for all subjects}

\item{margin_size}{A numeric vector of length *tens_dim* with the margin sizes of the response tensor}

\item{CNR}{The contrast-to-noise ratio, as defined by Welvaert, M. and Rosseel, Y., 2013}

\item{obs.var}{A value for the observational variance}
}
\value{
A list of length two with elements, *dat* and *truth*. The 'dat' object is a list with two elements, *Y* and *x*. *Y* is the tensor response in which the last two array margins represent the time step and the subject number, respectively. *x* is a single vector of length *n.time* that applies as the covariate in the tensor response regression. *truth* is the true coefficient tensor used in the data generation.
}
\description{
Bayesian Tensor Response Regression Data Generation
}
