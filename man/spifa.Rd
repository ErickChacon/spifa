% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spifa-package.R, R/spifa.R
\docType{package}
\name{spifa}
\alias{spifa-package}
\alias{spifa}
\title{spifa: Spatial (Geostatistical) Item Factor Analysis}
\usage{
spifa(
  responses,
  pred_formula = NULL,
  coords = NULL,
  data = NULL,
  nfactors,
  ngp = nfactors,
  niter = 1000,
  thin = 10,
  standardize = TRUE,
  constrains = list(A = NULL, W = NULL, V_sd = rep(1, nfactors)),
  adaptive = list(Sigma = NULL, Sigma_R = NULL, Sigma_gp_sd = NULL, Sigma_gp_phi = NULL,
    scale = 1, C = 0.7, alpha = 0.8, accep_prob = 0.234),
  c_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
  A_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
  R_opt = list(initial = NULL, prior_eta = 1.5),
  B_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
  sigmas_gp_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
  phi_gp_opt = list(initial = NULL, prior_mean = NULL, prior_sd = NULL),
  execute = TRUE
)
}
\arguments{
\item{constrains}{Named list of constrains associated to the factor model. Accepted
names are `A`, `W`, and `V_sd`. The restrictions in the discrimination paramater should
be place in the element `A` with same dimensions as the discrimination matrix (nitems x
nfactors). A value of 0 indicates that link betwen the item and the factor is disabled
and 1 indicates that it remain active and the coefficient associated will be estimated.
The restrictions for the gaussian processes should be place in the element `W` with
dimensions nfactors x ngps, such as a value of 0 indicates a link disconnected between
the factor and the gp while 1 indicates that it remains active. The restrictions with
respect to the standard deviation of the error terms inside the latent factors should
be placed in the element `V_sd`, which should be a vector (length nfactors) providing
the fixed values for the error standard deviation. If the model includes predictors or
gaussian processes, it is recomended to be lower than 1.}

\item{adaptive}{Named list of hyperparameters associated with the adaptive sampling.
The adaptive sampling is done jointly for the `correlation` parameters, `standard
deviation of the gps` and `scale parameter of the gps`. The matrix `Sigma` can be
provided as the full covariance matrix of these parameters for the proposal
distribution. Otherwise, part of this matrix can be provided by using the elements
`Sigma`, `Sigma_R`, `Sigma_gp_sd` and `Sigma_gp_phi`. Additional elements are `scale`,
`C`, `alpha` and `accep_prob` which are hyperparameters of the adaptive sampling
proposed in Andrieu and Thomas (2008).}

\item{A_opt}{Same as C_opt but for the discrimination parameters.}

\item{R_opt}{Same as C_opt but for the correlation parameters. This list only accepts
`initial` value and `prior_eta` associated to the LKJ prior.}

\item{B_opt}{Same as C_opt but for the regression parameters.}

\item{sigmas_gp_opt}{Same as C_opt but for the standard deviation of the gaussian
processes.}

\item{phi_gp_opt}{Same as C_opt but for the scale parameter of the gaussian processes.}

\item{execute}{Logical value to run sampler or not. TRUE by default.}

\item{C_opt}{Named list of initial values and hyperparameters for the easiness
parameters. The initial value is provided in the element `initial`, and the prior meand
and standard deviation are provided in the elements `prior_mean` and `prior_sd`
respectively.}
}
\value{
return.
}
\description{
The \code{spifa} package

\code{function} description.
}
\details{
details.
}
\examples{

bla

}
\author{
Erick A. Chacon-Montalvan \email{e.chaconmontalvan@lancaster.ac.uk}

Emanuele Giorgi, \email{e.giorgi@lancaster.ac.uk}

Benjamin M. Taylor, \email{b.taylor1@lancaster.ac.uk}

Erick A. Chacón-Montalván
}
