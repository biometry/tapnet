#' Fit the tapnet model to a network
#'
#' Estimates the parameters of the tapnet model by log-likelihood based on the oberved network(s)
#'
#' The core function for using the tapnet approach: it fits the model to the data (= networks). Then, the estimated parameters can be used to predict to other networks (using \code{\link{predict_tapnet}}).
#'
#' @aliases fit_tapnet
#' 
#' @param tapnet a tapnet object;
#' @param ini initial parameter values for the optimization; optional;
#' @param tmatch_type_pem type of trait matching function for latent traits, currently "normal" or "shiftlnorm";
#' @param tmatch_type_obs type of trait matching function for observed traits, currently "normal" or "shiftlnorm";
#' @param TmatchMatrixList list of independent trait-matching matrices (one per network);
#' @param lambda LASSO shrinkage factor for latent trait parameters;
#' @param method Optimization method (most derivative-based approaches will not work! SANN is a (slow) alternative to the default);
#' @param maxit Maximum number of steps for optimization;
#' @param hessian logical: output hessian for calculation of standard errors?
#' @param obj_function Objective function for the optimization, either "multinom" or "sq_diff" (or "bjorn");
#' @param fit.delta logical; should the parameter delta be fitted? It allows tapnet to down-weigh the importance of trait matching relative to abundances. Defaults to FALSE.
#' 
#' @return A tapnet-fit object, containing the tapnet model parameters as entries "par_opt", the settings of the tmatch_type for PEMs and observed traits, the parameter set for lambda, the optimisation method set, along with its maxit-value, and, finally, the output of the call to \code{optim}, including the target value (the negative log-likelihood), the convergence report and the parameters as fitted \emph{at the transformed scale}. Note that the entries under "opt" will not be the same as those under "par_opt"!
#' 
#' 
#' 
#' @references Benadi et al. in prep
#'
#' @author Gita Benadi <gita.benadi@biom.uni-freiburg.de> and Carsten Dormann <carsten.dormann@biom.uni-freiburg.de>
#'
#' @examples
#' \donttest{ # takes about 35 s
#'  data(Tinoco)
#'  tap <- make_tapnet(tree_low = plant_tree, tree_high = humm_tree, networks = networks[2:3], 
#'         traits_low = plant_traits, traits_high = humm_traits, npems_lat = 4)
#'  fit <- fit_tapnet(tap) # fits to networks 2 and 3 only
#'  str(fit)
#' }  
#'
#' @export
fit_tapnet <- function(tapnet, # a tapnet object
                       ini = NULL, # initial parameter values for the optimization
                       tmatch_type_pem = "normal", # type of trait matching function for latent traits
                       tmatch_type_obs = "normal", # type of trait matching function for observed traits
                       TmatchMatrixList = NULL, # list of pre-formulated trait-matching probability matrices, e.g. for phenology or known forbidden links; will be multiplied onto internally computed trait matching; must have as many matrices as there are networks in the tapnet object!
                       lambda = 0, # LASSO shrinkage factor for latent trait parameters
                       method = "Nelder", # Optimization method
                       maxit = 50000, # Maximum number of steps for optimization
                       hessian = FALSE, # Output hessian for calculation of standard errors?
                       obj_function = "multinom", # Objective function for the optimization,
                       # either "multinom" or "sq_diff"
                       fit.delta=FALSE # shall the abundance-modifier delta be fit?
) {
  # Fit a model of interaction probabilities based on traits, abundances and phylogenies
  # to one or more observed interaction networks
  
  if (!inherits(tapnet, "tapnet")) {
    stop("This is no tapnet object! Please use function make_tapnet or simulate_tapnet to create one.")
  }
  
  if (!is.null(TmatchMatrixList)) {
    if (length(TmatchMatrixList) != length(tapnet$networks) ) stop("TmatchMatrixList must have as many matrices as there are networks in the tapnet object!")
  }
  
  
  # if (tapnet$trees$low) ... # check whether a PEM-fit is asked for but no phylogeny is provided; if so, ask user to set tmatch_type_pem to "no" to ignore this information
  
  # Check length of initial values vector, create the vector if it hasn't been supplied
  pems_low <- lapply(tapnet$networks, function(x) x$pems[[1]])
  pems_high <- lapply(tapnet$networks, function(x) x$pems[[2]])
  pem_names_low <- unique(unlist(lapply(pems_low, colnames)))
  if (!is.null(pem_names_low)) pem_names_low <- pem_names_low[order(nchar(pem_names_low), pem_names_low)]
  pem_names_high <- unique(unlist(lapply(pems_high, colnames)))
  if (!is.null(pem_names_high)) pem_names_high <- pem_names_high[order(nchar(pem_names_high), pem_names_high)]
  
  if (is.null(tapnet$traits_all$low)) ntraits <- 0 else ntraits <- ncol(tapnet$traits_all$low)
  
  if (tmatch_type_pem =="no"){
    nparams <- 1 + ntraits + 1 # number of parameters to be fitted; one less without PEMs (the pem-shift parameter)
  } else {
    nparams <- length(pem_names_low) + length(pem_names_high) + 2 + ntraits + 1 # number of parameters to be fitted:
    # linear combination parameters for lower and higher level (as many as PEMs), one PEM shift parameter,
    # one PEM trait matching parameter, as many additional trait matching parameters as observed traits,
    # one abundance-weighting parameter (delta)
  }
  if (is.null(ini)) { # define here the default starting values for the optimiser!
    if (tmatch_type_obs == "normal") ini <- c(0, rep(1, nparams-1) ) # runif(nparams) # 1st parameter will be exponentiated
    if (tmatch_type_obs == "shiftlnorm") ini <- c(.1, rep(0.1, nparams-1) ) # runif(nparams)
    ini[length(ini)] <- 0 # approx. 0.5 = plogis(0) as value for delta
  } else {
    if (length(ini) != nparams) stop("Number of initial values must equal number of parameters to be fitted!") # with delta!
  }
  
  if (tmatch_type_pem == "no"){
    if (is.null(tapnet$traits_all$low)) {
      names(ini) <- c("tmatch_width_pem", "delta")
    } else {
      names(ini) <- c("tmatch_width_pem", paste0("tmatch_width_obs", 1:ncol(tapnet$traits_all$low)), "delta")
    }  
  } else {
    if (is.null(tapnet$traits_all$low)) {
      names(ini) <- c(pem_names_low, pem_names_high, "pem_shift", "tmatch_width_pem", "delta")
    } else {
      names(ini) <- c(pem_names_low, pem_names_high, "pem_shift", "tmatch_width_pem",
                    paste0("tmatch_width_obs", 1:ncol(tapnet$traits_all$low)), "delta")
    }
  }
  if (!fit.delta) ini <- ini[-which(names(ini) == "delta")]
  #test:
  # loglik_tapnet(ini, networks = tapnet$networks, tmatch_type_pem = tmatch_type_pem, tmatch_type_obs = tmatch_type_obs, lambda = lambda, fit.delta=F)
  
  # Optimization
  opt <- optim(par = ini, fn = loglik_tapnet, networks = tapnet$networks,
               tmatch_type_pem = tmatch_type_pem, tmatch_type_obs = tmatch_type_obs, lambda = lambda, fit.delta=fit.delta,
               control = list(maxit = maxit), method = method, hessian = hessian, obj_function = obj_function, TmatchMatrixList=TmatchMatrixList)

  # Convert optimized parameter vector to a named list 
  par_opt <- param_vec2list(opt$par, n = length(pem_names_low), m = length(pem_names_high), fit.delta=fit.delta)
  if (!fit.delta) par_opt[["delta"]] <- c("delta constant"=1)
  
  out <- list(par_opt = par_opt, tmatch_type_pem = tmatch_type_pem, tmatch_type_obs = tmatch_type_obs, lambda = lambda, method = method, maxit = maxit, opt = opt)
  class(out) <- "fitted.tapnet"
  attr(out, "tapnet_name") <- as.character(substitute(tapnet))
  return(out)
}