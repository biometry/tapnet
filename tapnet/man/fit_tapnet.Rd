% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_tapnet.R
\name{fit_tapnet}
\alias{fit_tapnet}
\title{Fit the tapnet model to a network}
\usage{
fit_tapnet(
  tapnet,
  ini = NULL,
  tmatch_type_pem = "normal",
  tmatch_type_obs = "normal",
  TmatchMatrixList = NULL,
  lambda = 0,
  method = "Nelder",
  maxit = 50000,
  hessian = FALSE,
  obj_function = "multinom",
  fit.delta = FALSE
)
}
\arguments{
\item{tapnet}{a tapnet object;}

\item{ini}{initial parameter values for the optimization; optional;}

\item{tmatch_type_pem}{type of trait matching function for latent traits, currently "normal" or "shiftlnorm";}

\item{tmatch_type_obs}{type of trait matching function for observed traits, currently "normal" or "shiftlnorm";}

\item{TmatchMatrixList}{list of independent trait-matching matrices (one per network);}

\item{lambda}{LASSO shrinkage factor for latent trait parameters;}

\item{method}{Optimization method (most derivative-based approaches will not work! SANN is a (slow) alternative to the default);}

\item{maxit}{Maximum number of steps for optimization;}

\item{hessian}{logical: output hessian for calculation of standard errors?}

\item{obj_function}{Objective function for the optimization, either "multinom" or "sq_diff" (or "bjorn");}

\item{fit.delta}{logical; should the parameter delta be fitted? It allows tapnet to down-weigh the importance of trait matching relative to abundances. Defaults to FALSE.}
}
\value{
A tapnet-fit object, containing the tapnet model parameters as entries "par_opt", the settings of the tmatch_type for PEMs and observed traits, the parameter set for lambda, the optimisation method set, along with its maxit-value, and, finally, the output of the call to \code{optim}, including the target value (the negative log-likelihood), the convergence report and the parameters as fitted \emph{at the transformed scale}. Note that the entries under "opt" will not be the same as those under "par_opt"!
}
\description{
Estimates the parameters of the tapnet model by log-likelihood based on the oberved network(s)
}
\details{
The core function for using the tapnet approach: it fits the model to the data (= networks). Then, the estimated parameters can be used to predict to other networks (using \code{\link{predict_tapnet}}).
}
\examples{
\donttest{ # takes about 35 s
 data(Tinoco)
 tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], 
        traits_low = plant_traits, traits_high = humm_traits, npems_lat = 4)
 fit <- fit_tapnet(tap) # fits to networks 2 and 3 only
 str(fit)
}  

}
\references{
Benadi et al. in prep
}
\author{
Gita Benadi <gita.benadi@biom.uni-freiburg.de> and Carsten Dormann <carsten.dormann@biom.uni-freiburg.de>
}
